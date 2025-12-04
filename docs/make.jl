using Documenter, UsefulFunctions

makedocs(
    modules=[UsefulFunctions],
    sitename="Useful Functions",
    pages=Any["Home"=>"index.md","Indices"=>"indices.md",hide("Root Finding"=>"root_finding.md")]
    )

deploydocs(
    repo="github.com/biplab37/UsefulFunctions.jl.git"
)

