module Input

using JuES

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
    mol = psi4.geometry(cont["molecule"])
    psi4_options = Dict{Any,Any}(
                       "scf_type"=>"pk"
    )
    JuES_options = Dict{Any,Any}(
                       :doprint=>false,
                       :maxit=>40,
                       :return_T=>false,
                       :quiet=>true
                      )
    for option in keys(cont["options"]) #sort through options to discern psi4 keywords vs JuES keywords
        if option in psi4_opt_list
            psi4_options[option] = cont["options"][option]
        else
            JuES_options[option] = cont["options"][option]
        end
    end
    if JuES_options[:quiet] == true
        psi4.core.be_quiet()
    end
    psi4.set_options(psi4_options)
    e,wfn = psi4.energy("scf",return_wfn=true)
    JuWfn = JuES.Wavefunction.Wfn(wfn)
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
