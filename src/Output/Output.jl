module Output
import JuES
using Printf
using Formatting
import Printf.decode_dec
import Printf.fix_dec
import Printf.print_fixed
#using JuES.Options

#printstyle = JuES.Options.printstyle
printstyle = "stdout"
export output
export @output

if printstyle == "stdout"
    function output(str,x...)
        f = format(str,x...)
        print(f)
    end
    macro output(str,x...)
        return quote
            local f = format($str, $(x...))
            print(f)
        end
    end
elseif printstyle == "file"
    function output(str,x...)
        f = format(str,x...)
        open("output.dat","a") do file
            write(file,f)
            flush(file)
        end
    end
    macro output(str,x...)
        return quote
            local f = format($str, $(x...))
            open("output.dat","a") do file
                write(file,f)
                flush(file)
            end
        end
    end
elseif printstyle == "both"
    function output(str,x...)
        f = format(str,x...)
        print(f)
        open("output.dat","a") do file
            write(file,f)
            flush(file)
        end
    end
    macro output(str,x...)
        if length(x) >= 1
            return quote
                local f = format($str, $([esc(i) for i in x]...))
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
                print(f)
            end
        else
            return quote
                local f = format($str)
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
                print(f)
            end
        end
    end
elseif printstyle == "none"
    function output(str,x...)
    end
    macro output(str,x...)
    end
end

end #module
