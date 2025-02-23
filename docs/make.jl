using StochasticDominance
using Documenter

makedocs(
	sitename =  "StochasticDominance.jl",
	authors = "Rajmadan Lakshmanan",
	clean = true,
	doctest = false,
	format = Documenter.HTML(),
	pages = ["Home" => "index.md",
		"Tutorials" => Any["tutorial/tutorial1.md",
				    "tutorial/tutorial2.md",
					"tutorial/tutorial3.md",
					"tutorial/tutorial4.md",
				    ],
		"API Reference" => "api.md"
		]
)

deploydocs(
    repo = "github.com/rajmadan96/StochasticDominance.jl.git",
    branch = "gh-pages",
    push_preview = true,
    forcepush = true,
    deploy_config = Documenter.GitHubActions(),  # Ensure it's using GitHub Actions
)





