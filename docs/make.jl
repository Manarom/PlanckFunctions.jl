include(joinpath(@__DIR__, "..", "test", "tests data", "TestingData.jl"))

using Documenter, PlanckFunctions, .TestingData

mathengine = Documenter.MathJax3()

makedocs(
    sitename = "PlanckFunctions.jl",
    highlightsig = false,
    checkdocs = :none,
    format = Documenter.HTML(
        size_threshold = 2000 * 2^10,
        mathengine = mathengine
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Examples" => "pluto_tests_git.md"
        ],
        "Modules" => [
            "PlanckFunctions" => "PlanckFunctions.md",
            "TestingData" => "TestingData.md"
        ] 
    ]
)

deploydocs(;
    devbranch = "master",
    repo = "://github.com",
    devurl = "dev",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"]
)