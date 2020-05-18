module Output
import JuES
using Formatting
#using JuES.Options

#printstyle = JuES.Options.printstyle
if !isdefined(JuES.Output,:printstyle)
    printstyle = ["stdout"]
end
export output
export @output
function set_print(pstyle)
    if pstyle in ["none","file","stdout","both"]
        JuES.Output.printstyle[1] = pstyle
        #revise(JuES.Output)
    else
        return false
    end
end

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
    return quote
        if JuES.Output.printstyle[1] == "stdout"
                local f = format($str, $([esc(i) for i in x]...))
                print(f)
        elseif JuES.Output.printstyle[1] == "file"
                local f = format($str, $([esc(i) for i in x]...))
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
        elseif JuES.Output.printstyle[1] == "both"
            if length(x) >= 1
                local f = format($str, $([esc(i) for i in x]...))
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
                print(f)
            else
                local f = format($str)
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
            end
        elseif JuES.Output.printstyle[1] == "none"
        end
    end
end
end #module
