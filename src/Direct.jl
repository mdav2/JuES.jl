module Direct
export ao
export func_shell
using PyCall
using Wavefunction
const psi4 = PyNULL()
function __init__()
	copy!(psi4,pyimport("psi4"))
end

function ao_function(wfn::DirectWfn,p,q,r,s)
	"""
	Computes a single ERI in AO basis
	Very wasteful (throws out rest of shell)
	use for debugging
	"""
	shell = ao_shell(wfn,p,q,r,s)
	return shell[po,qo,ro,so]
end
function func_shell(wfn::PyObject)
	"""
	Matches basis function number and corresponding shell number
	"""
	basis = wfn.basis
	nao = basis.nao()
	output = Dict()
	for i in 1:1:nao
		x = basis.function_to_shell(i-1)
		if haskey(output, x)
			push!(output[x],i)
		else
			output[x] = [i]
		end
	end
	return output
end
function ao_shell(wfn::DirectWfn,p,q,r,s)
	"""
	Computes the AO eri shell (pq|rs)
	"""
	basis = wfn.basis
	mints = wfn.mints
	ps = basis.function_to_shell(p-1)
	qs = basis.function_to_shell(q-1)
	rs = basis.function_to_shell(r-1)
	ss = basis.function_to_shell(s-1)

	po = p - basis.shell_to_ao_function(ps)# + 1
	qo = q - basis.shell_to_ao_function(qs)# + 1
	ro = r - basis.shell_to_ao_function(rs)# + 1
	so = s - basis.shell_to_ao_function(ss)# + 1
	shell = mints.ao_eri_shell(ps,qs,rs,ss).to_array()
	return shell
end
function direct_MP2(wfn::DirectWfn)
	"""
	Algorithm
	--
	for each occ <ij|-->, compute partially transformed
	integrals <ij|λσ>, store in memory

	for each vir <--|ab> contract <ij|λσ> ⋅ C[a,λ] ⋅ C[b,σ]
	and its transpose <--|ba>
	
	compute contribution of <ij|ab> to δE(MP2)

	   12 34      12 34     12 34        
	= <ij|ab> ( 2<ij|ab> - <ij|ba> ) / Δ[ijab]
	= <ij|ab> ( 2<ij|ab> - <ij|ba> ) / Δ[ijab]
	------------------------------------------
	   13 24      13 24     13 24        
	= [ia|jb] ( 2[ia|jb] - [ib|ja] ) / Δ[ijab]
	"""
	nshell = basis.nshell()
	#figure out how to make use of shell/batches
	# "         how to make partially transformed integrals
	# "         how to do contraction
	#for storing partially transformed integrals
	partial = zeros(Float64,wfn.nalpha,wfn.nalpha,wfn.nmo,wfn.nmo)
    for s in R
        for si in R
            for rh in R
                for nu in R
                    for mu in R
                        g1[mu,nu,rh,s] += gao[mu,nu,rh,si]*C[si,s]
                    end
                end
            end
        end
    end
    for s in R
        for r in R
            for rh in R
                for nu in R
                    for mu in R
                        g2[mu,nu,r,s] += g1[mu,nu,rh,s]*C[rh,r]
                    end
                end
            end
        end
    end
end
end
