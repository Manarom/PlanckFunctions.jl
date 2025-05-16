push!(LOAD_PATH,"../src/")
include("../test/tests data/TestingData.jl")
using Documenter,PlanckFunctions,.TestingData
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "PlanckFunctions.jl",
        repo="https://github.com/Manarom/PlanckFunctions.jl/blob/{commit}{path}#{line}",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = 2000 * 2^10),
        pages=[
                "Home" => "index.md"
                "Examples"=>["Examples" =>"pluto_tests_git.md"
                ]
                "Modules" => [
                    "PlanckFunctions" =>"PlanckFunctions.md"
                    "TestingData"=>"TestingData.md"
                 ] 
               ]#
			   )
deploydocs(;
                repo="github.com/Manarom/PlanckFunctions.jl", 
                devbranch = "master",
                devurl="dev",
                target = "build",
                branch = "gh-pages",
                versions = ["stable" => "v^", "v#.#" ]
        )