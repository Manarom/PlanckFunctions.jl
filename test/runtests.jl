import  PlanckFunctions as PF
using Test #,LinearAlgebra,Statistics
using ForwardDiff
include(joinpath(@__DIR__,"tests data","TestingData.jl")) # TestingData.benchmark_data
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
dfdx(f,x) = ForwardDiff.derivative(f,x)
d2fdx2(f,x) = dfdx(t->dfdx(f,t) , x)

@testset "PlanckFunctions.jl" begin
    # Write your tests here.
    for point in values(TestingData.benchmark_data) #iterating over data for various temperatures
        @show T = point.temperature
        @test T ≈ PF.temperature(PF.ibb(1.0 , T),1.0)

       print("Planck function tests ... ")

            _amat = fill(0.0 , (10,3))
            _lam = range(1.0,2.0,10)
            PF.a₁₂₃!(_amat , _lam , T)
            @test PF.ibb.(_lam , T) ≈ PF.ibb(_lam , _amat)
            _ibb = similar(collect(_lam))
            PF.ibb!(_ibb , _lam , T)
            @test _ibb ≈ PF.ibb.(_lam , T)
            _i_mat =  PF.ibb(_lam , [1000.0 , T])
            @test _i_mat[:,end] ≈ _ibb
            _ibb = similar(collect(_lam))
            PF.ibb!(_ibb , _lam , _amat)
            @test PF.ibb.(_lam , T) ≈ _ibb

        println("ok")

        print("planck function temperature and lambda derivatives vs ForwardDiff ... ") 

            f_t = t->PF.ibb(2.0 , t)
            f_l = l->PF.ibb(l , T)
            @test PF.∇ₜibb(2.0 , T) ≈ dfdx(f_t , T)
            @test PF.∇²ₜibb(2.0 , T) ≈ d2fdx2(f_t , T)
            @test PF.∇ₗibb(2.0 , T) ≈ dfdx(f_l , 2.0)
            @test PF.∇²ₗibb(2.0 , T) ≈ d2fdx2(f_l , 2.0)

        println("ok")

        print("various in-place versions ... ") 

            # ∇ₜibb!(g::AbstractMatrix , λ::AbstractVector , T::AbstractVector)
            _i_mat =  fill(0.0 , (length(_lam) , 3))
            PF.∇ₜibb!(_i_mat , _lam , [500.0 , 1000.0 , T])
            @test _i_mat ≈ PF.∇ₜibb(_lam , [500.0 , 1000.0 , T])
            # ∇ₜibb!(g::AbstractVector , λ::AbstractVector , T , amat::AbstractMatrix)
            _ibb = similar(collect(_lam))
            _amat = fill(0.0 , (length(_lam),3))
            PF.a₁₂₃!(_amat , _lam , T)
            PF.∇ₜibb!(_ibb , _lam , T , _amat)
            @test _ibb ≈ PF.∇ₜibb.(_lam , T)

            _i_mat =  fill(0.0 , (length(_lam) , 3))
            PF.∇²ₜibb!(_i_mat , _lam , [500.0 , 1000.0 , T])
            @test _i_mat ≈ PF.∇²ₜibb(_lam , [500.0 , 1000.0 , T])
            # ∇²ₜibb!(h::AbstractVector , λ::AbstractVector , T::Number , amat::AbstractMatrix)
            _ibb = similar(collect(_lam))
            _amat = fill(0.0 , (length(_lam),3))
            PF.a₁₂₃!(_amat , _lam , T)
            PF.∇²ₜibb!(_ibb , _lam , T,  _amat)
            @test _ibb ≈ PF.∇²ₜibb.(_lam , T)

        println("ok")

        print("combined functions ... ")

            @test all((PF.ibb(2.0 , T) , PF.∇ₗibb(2.0 , T) , PF.∇²ₗibb(2.0 , T))  .≈  PF.Dₗibb(2.0 , T))
            _lam = range(1.0,2.0,10)
            _T = [1234.6 , 456.7 , T]
            @test all((PF.ibb(_lam , _T) , PF.∇ₜibb(_lam , _T) , PF.∇²ₜibb(_lam , _T))  .≈  PF.Dₜibb(_lam , _T))
            _i = ntuple(_->fill(0.0 ,length(_lam)) , 3)
            @test all((PF.ibb.(_lam , T) , PF.∇ₜibb.(_lam , T) , PF.∇²ₜibb.(_lam , _T))  .≈  PF.Dₜibb(_lam , T))

        println("ok")

        print("Planck function power and band integration and derivatives with temperature ... ")

            @test PF.∇ₜpower(T) ≈ dfdx(PF.power , T)
            @test PF.∇²ₜpower(T) ≈ d2fdx2(PF.power , T)
            bp(T) = PF.band_power(T , λₗ = 1.0 , λᵣ=6.0 , tol =1e-12)
            @test PF.∇ₜband_power(T , λₗ = 1.0 , λᵣ=6.0 ) ≈ dfdx(bp , T)
            @test PF.∇²ₜband_power(T , λₗ = 1.0 , λᵣ=6.0 ) ≈ d2fdx2(bp , T) atol=1e-6 # there is a difference wiht ForwardDiff about 1e-8
            # f_pow = PF.power()

        println("ok")

        print("testing based on externally provided data ... ") 

            λₗ = point.lower #lower wavelength
            λᵣ = point.upper #upper wavelength
            sp = point.spectrum 
            # total radiance test
            @test PF.power(T) ≈ point.total_radiance rtol=1e-4 # (possible due to some difference in Stefan constant there is a discrepancy)
            #peak wavelength test
            @test PF.λₘ(T) ≈ point.peak_wavelength rtol=1e-6 
            #simple Planck function
            @test PF.ibb(PF.λₘ(T),T) ≈ point.peak_spectral_radiance rtol=1e-6
            #testing the in-band radiance
            @test PF.band_power(T,λₗ =λₗ ,λᵣ =λᵣ ) ≈ point.band_radiance rtol=1e-4
            # testing least-square difference between calculated spectra
            points_number = size(sp,1)
            n = sqrt(sum(i->i^2, PF.ibb.(sp[:,1],T) .-sp[:,2]))/points_number
            @test n ≈ 0 atol=1e-2

        println("ok")
    end
end
