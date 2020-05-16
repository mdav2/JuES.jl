module Input

using JuES
using Serialization
using PyCall: PyObject

export run #shorthand for read->exec
export read
export exec

include("psi4keywords.dat")

"""
    run(fname::String)

read and execute the input file specified by fname. 
"""
function run(fname::String = "input.dat")
    cont = read(fname)
    exec(cont)
end

"""
    read(fname::String)

read an input file and convert to a Dict. exec takes the output of this function as input.
"""
function read(fname)
    if isfile(fname)
        lines = readlines(fname)
    else
        return false
    end
    molecule,options,command = parse(lines)
    cont = Dict("molecule"=>molecule,"options"=>options,"command"=>command)
end

"""
    exec(cont::Dict)

execute the computation specified in cont.
cont is composed of:
    Dict("molecule"=>molecule::String,
         "options"=>options::Dict,
         "command"=>command::function)
"""
function exec(cont)
    if !occursin("symmetry c1",lowercase(cont["molecule"]))
        cont["molecule"] *= "symmetry c1"
    end
    psi4_options = Dict{Any,Any}(
                       "scf_type"=>"pk"
    )
    JuES_options = Dict{Any,Any}(
                       :doprint=>false,
                       :maxit=>40,
                       :return_T=>false,
                       :quiet=>true,
                       :T=>Float64
                      )
    for option in keys(cont["options"]) #sort through options to discern psi4 keywords vs JuES keywords
        if option in psi4_opt_list
            psi4_options[option] = cont["options"][option]
        else
            JuES_options[option] = cont["options"][option]
        end
    end
    T = JuES_options[:T]
    PATH_TO_PSI4 = joinpath(dirname(pathof(JuES)),"Psi4/Psi4.jl")
    PATH_TO_INPUT = "/tmp/input.jl"
    PATH_TO_SCR = "/tmp/wfn.jld"
    makey = """
    using Serialization
    #using JLD
    #using JuES.Wavefunction: Wfn
    include("$PATH_TO_PSI4")
    mol = Psi4.psi4.geometry(\"\"\"$(cont["molecule"])\"\"\")
    Psi4.psi4.core.be_quiet()
    Psi4.psi4.set_options($psi4_options)
    e,wfn = Psi4.psi4.energy("scf",return_wfn=true)
    mints = Psi4.psi4.core.MintsHelper(wfn.basisset())
    dt = $T
    energy = wfn.energy()
    dummy2 = Array{dt}(undef, 0, 0) #placeholder 2D array
    dummy4 = Array{dt}(undef, 0, 0, 0, 0) #placeholder 4D array
    vnuc = wfn.molecule().nuclear_repulsion_energy()
    basis = wfn.basisset()
    nbf = wfn.nmo()
    nocca = wfn.nalpha()
    nvira = nbf - nocca
    noccb = wfn.nbeta()
    nvirb = nbf - noccb
    epsa = convert(Array{dt,1}, wfn.epsilon_a().to_array()) #orbital eigenvalues
    epsb = convert(Array{dt,1}, wfn.epsilon_b().to_array()) #orbital eigenvalues
    _Ca = wfn.Ca() #as psi4 Matrix objects for use with MintsHelper
    _Cb = wfn.Cb() #as psi4 Matrix objects for use with MintsHelper
    _Cao = wfn.Ca_subset("AO", "OCC")
    _Cav = wfn.Ca_subset("AO", "VIR")
    _Cbo = wfn.Cb_subset("AO", "OCC")
    _Cbv = wfn.Cb_subset("AO", "VIR")
    Cao = convert(Array{dt,2}, wfn.Ca_subset("AO", "OCC").to_array())
    Cav = convert(Array{dt,2}, wfn.Ca_subset("AO", "VIR").to_array())
    Cbo = convert(Array{dt,2}, wfn.Cb_subset("AO", "OCC").to_array())
    Cbv = convert(Array{dt,2}, wfn.Cb_subset("AO", "VIR").to_array())
    Ca = convert(Array{dt,2}, _Ca.to_array())
    Cb = convert(Array{dt,2}, _Cb.to_array())
    hao = convert(Array{dt,2}, wfn.H().to_array()) #core hamiltonian in AO
    ao_eri = convert(Array{dt,4}, mints.ao_eri().to_array()) #AO basis integrals

    
    odict = Dict("energy" => energy,
                 "vnuc" => vnuc,
                 "nocca" => nocca,
                 "noccb" => noccb,
                 "nvira" => nvira,
                 "nvirb" => nvirb,
                 "nbf" => nbf,
                 "unrestricted" => false,
                 "Ca" => Ca,
                 "Cb" => Cb,
                 "Cao" => Cao,
                 "Cbo" => Cbo,
                 "Cav" => Cav,
                 "Cbv" => Cbv,
                 "hao" => hao,
                 "epsa" => epsa,
                 "epsb" => epsb,
                 "ao_eri" => ao_eri,
                 "T" => dt)

    serialize("/tmp/jues.1",odict)
    #save("$PATH_TO_SCR","wfn",wfn)
    """
    open(PATH_TO_INPUT,"w") do f
        write(f,makey)
    end
    Base.run(`julia $PATH_TO_INPUT`)
    wfn = deserialize("/tmp/jues.1")
    JuWfn = JuES.Wavefunction.Wfn{Float64}(wfn)
    com = cont["command"]
    E = com(JuWfn; JuES_options...)
    return E
end

function parse(lines;delim="--")
    chunks = split(join(lines,"\n"),delim)
    for chunk in chunks
        expr = Meta.parse(chunk)
        eval(expr)
    end
    molecule,options,command
end


end # module
