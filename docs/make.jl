using StochasticDominance
using Documenter

makedocs(
	sitename =  "StochasticDominance.jl",
	authors = "Rajmadan Lakshmanan",
	clean = true,
	doctest = false,
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	pages = ["Home" => "index.md",
		"Tutorials" => Any["tutorial/tutorial1.md",
				    "tutorial/tutorial2.md",
				    ]
		]
)

deploydocs(; 
	repo = "git@github.com:rajmadan96/StochasticDominance.jl.git",
	target="build",
	#devbranch = "main",
    	push_preview=false
)


