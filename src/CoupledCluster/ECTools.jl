module ECTools
using JuES
using JuES.Output

export read_ci_dets

function read_ci_dets(fname::String; max_exc::Int=4)
    open(fname) do data
        for ln in eachline(data)
            m = match(r"\s*\*\d*\s+([-]?\d+\.\d+)\s+\(.+,.+\)\s+(.+)", ln)
        end
    end
end
    
end #Module
