module Output
import JuES
using Printf
using Formatting
import Printf.decode_dec
import Printf.fix_dec
import Printf.print_fixed
#using JuES.Options

#printstyle = JuES.Options.printstyle
printstyle = ["none"]
export output
export @output

function output(str,x...)
    if JuES.Output.printstyle[1] == "stdout"
        f = format(str,x...)
        print(f)
    elseif JuES.Output.printstyle[1] == "file"
        f = format(str,x...)
        open("output.dat","a") do file
            write(file,f)
            flush(file)
        end
    elseif JuES.Output.printstyle[1] == "both"
        f = format(str,x...)
        print(f)
        open("output.dat","a") do file
            write(file,f)
            flush(file)
        end
    elseif JuES.Output.printstyle[1] == "none"
    end
end
macro output(str,x...)
    if JuES.Output.printstyle[1] == "stdout"
        return quote
            local f = format($str, $([esc(i) for i in x]...))
            print(f)
        end
    elseif JuES.Output.printstyle[1] == "file"
        return quote
            local f = format($str, $([esc(i) for i in x]...))
            open("output.dat","a") do file
                write(file,f)
                flush(file)
            end
        end
    elseif JuES.Output.printstyle[1] == "both"
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
            end
        end
    elseif JuES.Output.printstyle[1] == "none"
    end
end
end #module
