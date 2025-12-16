using Documenter, UsefulFunctions

makedocs(
    modules=[UsefulFunctions],
    sitename="Useful Functions (by Biplab)",
    pages=["Home"=>"index.md",
        "Functions"=>
        [
            "Root Finding"=>"root_finding.md",
            "Differential Equations"=>"diffeqn.md",
            "Integrations"=>"integration.md",
            "Interpolations"=>"interpolation.md",
            "Filters"=>"filters.md",
            "Miscellaneous"=>"miscellaneous.md"
        ],
        "Indices"=>"indices.md"
        ]
    )

deploydocs(
    repo="github.com/biplab37/UsefulFunctions.jl.git"
)

