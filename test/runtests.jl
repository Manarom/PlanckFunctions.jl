using PlanckFunctions
using Test,LinearAlgebra,Statistics
include(joinpath(@__DIR__,"tests data","testing_data.jl")) # TestingData.benchmark_data
#=
temperature # Kelvins
lower # lower wavelength limit μm
upper # upper wavelength limit μm
total_radiance # W/(m²⋅sr)  
peak_wavelength # μm
peak_spectral_radiance # W/m²/sr/μm
band_radiance # W/m²/sr total radiance within the band lower-
spectrum # 
=#

@testset "PlanckFunctions.jl" begin
    # Write your tests here.
    for d in TestingData.benchmark_data #iterating over data for various temperatures
        point = d[2]
        @show T = point.temperature
        λₗ = point.lower #lower wavelength
        λᵣ = point.upper #upper wavelength
        sp = point.spectrum 
        @test PlanckFunctions.power(T) ≈ point.total_radiance rtol=1e-4 # total radiance (possible due to some difference in Stefan constant )
        @test PlanckFunctions.λₘ(T) ≈ point.peak_wavelength rtol=1e-6 # peak wavelength
        @test PlanckFunctions.ibb(PlanckFunctions.λₘ(T),T) ≈ point.peak_spectral_radiance rtol=1e-6
        @test band_power(T,λₗ =λₗ ,λᵣ =λᵣ ) ≈ point.band_radiance rtol=1e-4
        @show points_number = size(sp,1)
        n = norm(ibb.(sp[:,1],T) .-sp[:,2])/points_number
        @test n ≈ 0 atol=1e-2
    end
end
