push!(LOAD_PATH,"../src/")

using ThreeBodyTB, Documenter
#using DocumenterLaTeX
    
current=pwd()
include("../deps/build.jl")
cd(current)

makedocs(sitename="ThreeBodyTB.jl Documentation")
#makedocs(sitename="ThreeBodyTB.jl Documentation", format = LaTeX() )


#push!(LOAD_PATH,"../src/")


#makedocs(sitename="ThreeBodyTB.jl Documentation")


#        prettyurls = get(ENV, "CI", nothing) == "true",
#        canonical = "https://oxfordcontrol.github.io/COSMO.jl/stable/",
#        assets = ["assets/favicon.ico"; "assets/github_buttons.js"; "assets/custom.css"],
#        analytics = "UA-134239283-1",
#  ),




@info "Making documentation..."
makedocs(
    sitename="ThreeBodyTB.jl Documentation",
    authors = "Kevin F. Garrity",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico", "assets/nist-combined.css"],
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => Any[
            "Running Calculations" =>  "ug_run.md",
            "Fit Coefficients" => "ug_fit.md",
        ],
        "Core User Interface" => Any[
            "Structs" => "structs.md",
            "Functions" => "core.md",
        ],
        "Compile/Python" => "compile.md",
        "Additional Docstrings" => "every.md"
    ]
)


@info "edit docs"


DD = ThreeBodyTB.DOCSDIR

stuff=readlines("$DD/nist_stuff/html_stuff.txt")
stuff2=readlines("$DD/nist_stuff/html_stuff_v2.txt")

function fix_html(f, s)
    lines=readlines(f)

    f2 = open(f, "w")
    for line in lines
        if occursin("</head>", line )
            line2 = replace(line, "</head>" => s[1]*"</head>")
            write(f2, line2*"\n")
        else
            write(f2, line*"\n")
        end

    end
    close(f2)

end
    
for d in readdir("$DD/build")
    if isdir("$DD/build/$d")
        f = "$DD/build/$d/index.html"
        if isfile(f)
            fix_html(f, stuff2)
        end
        if !isfile("$DD/build/$d/nist-combined.css")
            cp("$DD/nist_stuff/nist-combined.css", "$DD/build/$d/nist-combined.css")
        end
    end

end

fix_html("$DD/build/index.html", stuff)
if !isfile("$DD/build/nist-combined.css")
    cp("$DD/nist_stuff/nist-combined.css", "$DD/build/nist-combined.css")
end


@info "Deploy docs ..."

#deploydocs(;
#    repo="github.com/kfgarrity/ThreeBodyTB.jl",
#    push_preview=true,
#    devbranch = "main"
#)



#deploydocs()

#deploydocs(
#    repo = "github.com/kfgarrity/ThreeBodyTB.jl.git",
#    push_preview=true,
#)

