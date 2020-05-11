module Output
import JuES
using Printf
import Printf.decode_dec
import Printf.fix_dec
import Printf.print_fixed
#using JuES.Options

#printstyle = JuES.Options.printstyle
printstyle = "stdout"
export @output

if printstyle == "stdout"
    macro output(str,x...)
        if length(x) >= 1
            :(@printf $((str)) $(esc(x...)))
        else
            :(@printf $((str)) )#$(esc(x[1:end])))
        end

    end
elseif printstyle == "file"
    macro output(str,x...)
        if length(x) >= 1
            :(open("output.dat","a") do f
                  write(f,@sprintf $str $(esc(x...)))
              end)
        else
            :(open("output.dat","a") do f
                  write(f,@sprintf $str)
              end)
        end
    end
elseif printstyle == "both"
    macro output(str,x...)
        if length(x) >= 1
            :(open("output.dat","a") do f
                  write(f,@sprintf $str $(esc(x...)))
                  @printf $((str)) $(esc(x...))
              end)
        else
            :(open("output.dat","a") do f
                  write(f,@sprintf $str)
                  @printf $((str)) 
              end)
        end
    end
elseif printstyle == "none"
    macro output(str,x...)
        return nothing
    end
end

end #module
