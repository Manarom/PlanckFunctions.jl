push!(LOAD_PATH,"../src/")
include("../test/tests data/TestingData.jl")
using Documenter,PlanckFunctions,.TestingData
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "PlanckFunctions.jl",
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
