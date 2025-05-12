push!(LOAD_PATH,"../src/")
include("../src/PlanckFunctions.jl")
using Documenter,.PlanckFunctions
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "PlanckFunctions.jl",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = 2000 * 2^10),
        pages=[
                "Home" => "index.md"
                "Examples"=>["PlanckFunctions" =>"pluto_tests_git.md"
                ]
                "Modules" => [
                    "PlanckFunctions" =>"planck.md"
                 ] 
               ]#
			   )
