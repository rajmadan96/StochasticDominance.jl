using StochasticDominance
using Documenter
#push!(LOAD_PATH,"/home/mathe/lvhome10/users/personal/lraj/Nextcloud/SD/HOME1/users/personal/lraj/Nextcloud/SD/Github/StochasticDominance/Project.toml")
#DocMeta.setdocmeta!(StochasticDominance, :DocTestSetup, :(using StochasticDominance); recursive=true)

makedocs(
	sitename =  "StochasticDominance.jl",
	authors = "Rajmadan Lakshmanan",
	clean = true,
	doctest = false,
	format = Documenter.HTML(;canonical="https://rajmadan96.github.io/StochasticDominance.jl"),
	pages = ["Home" => "index.md",
		"Tutorials" => Any["tutorial/tutorial1.md",
				    "tutorial/tutorial2.md",
				    ]
		]
)

deploydocs(; 
	repo = "https://github.com/rajmadan96/StochasticDominance.jl.git",
	target = "build",
	devbranch = "master",
    	versions = ["stable" => "v^", "v#.#"]
)


