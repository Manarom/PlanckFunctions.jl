#debug_script

using Revise,BenchmarkTools

includet(joinpath(@__DIR__(),"PlanckFunctions.jl"))

l = collect(1.0:1e-2:30);
T = 2345.0

using .PlanckFunctions
#@which ibb(l,T)
@benchmark ibb.(l,T)
@benchmark ibb(l,T)

@benchmark ∇²ₗibb.(l,T)