import  PlanckFunctions as PF
using Test #,LinearAlgebra,Statistics
using ForwardDiff , QuadGK , ChainRulesCore
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

        print("planck function temperature and lambda derivatives and chain rules vs ForwardDiff ... ") 

            f_t = t->PF.ibb(2.0 , t)
            f_l = l->PF.ibb(l , T)
            r_t = t -> PF.spectral_ratio(0.8 , 1.1 , t)

            @test PF.∇ₜibb(2.0 , T) ≈ dfdx(f_t , T)
            @test PF.∇²ₜibb(2.0 , T) ≈ d2fdx2(f_t , T)
            @test PF.∇ₗibb(2.0 , T) ≈ dfdx(f_l , 2.0)
            @test PF.∇²ₗibb(2.0 , T) ≈ d2fdx2(f_l , 2.0)
            @test PF.∇ₜspectral_ratio(0.8 , 1.1 , T) ≈ dfdx(r_t , T)
            @test PF.∇²ₜspectral_ratio(0.8 , 1.1 , T) ≈ d2fdx2(r_t , T)
            (i1 , i2 , i3) = (r_t(T) , PF.∇ₜspectral_ratio(0.8 , 1.1 , T) , PF.∇²ₜspectral_ratio(0.8 , 1.1 , T))
            @test all( PF.Dₜspectral_ratio(0.8 , 1.1 , T) .≈ (i1 , i2 , i3))

            # function to check the derivatives of spectral band ratio's
            band_ration(T) = PF.band_power(T , λₗ = 1.0 , λᵣ=6.0 )/PF.band_power(T , λₗ = 4.0 , λᵣ=8.0 )
            @test PF.∇ₜspectral_band_ratio((1.0 , 6.0) , (4.0 , 8.0) , T) ≈ dfdx(band_ration , T) atol = 1e-5
            @test PF.∇²ₜspectral_band_ratio((1.0 , 6.0) , (4.0 , 8.0) , T) ≈ d2fdx2(band_ration , T) atol = 1e-5

            (i1 , i2 , i3) = (band_ration(T) , 
                            PF.∇ₜspectral_band_ratio((1.0 , 6.0) , (4.0 , 8.0) , T) , 
                            PF.∇²ₜspectral_band_ratio((1.0 , 6.0) , (4.0 , 8.0) , T))

            @test all( PF.Dₜspectral_band_ratio((1.0 , 6.0) , (4.0 , 8.0) , T) .≈ (i1 , i2 , i3))

            # testing forward differentiation chainrule for ibb 
            nt = ChainRulesCore.NoTangent()
            @test all((PF.ibb(2.0 , T) , PF.∇ₜibb(2.0 , T)) .≈ ChainRulesCore.frule(
                (nt , nt , 1.0), PF.ibb , 2.0 , T
            ))

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
            @test all((PF.ibb.(_lam , T) , PF.∇ₜibb.(_lam , T) , PF.∇²ₜibb.(_lam , T))  .≈  PF.Dₜibb(_lam , T))

        println("ok")

        print("Planck function power and band_power and derivatives with temperature ... ")

            @test PF.∇ₜpower(T) ≈ dfdx(PF.power , T)
            @test PF.∇²ₜpower(T) ≈ d2fdx2(PF.power , T)
            bp(T) = PF.band_power(T , λₗ = 1.0 , λᵣ=6.0 , tol =1e-12)
            @test PF.∇ₜband_power(T , λₗ = 1.0 , λᵣ=6.0 ) ≈ dfdx(bp , T)
            @test PF.∇²ₜband_power(T , λₗ = 1.0 , λᵣ=6.0 ) ≈ d2fdx2(bp , T) atol=1e-6 # there is a difference wiht ForwardDiff about 1e-8
            # f_pow = PF.power()
            bp_value = bp(T)
            @test PF.∇ₜband_power(T  , bp_value ; λₗ = 1.0 , λᵣ=6.0 ) ≈ dfdx(bp , T)
            @test PF.∇²ₜband_power(T  , bp_value ; λₗ = 1.0 , λᵣ=6.0 ) ≈ d2fdx2(bp , T) atol=1e-6

            

            @test all((
                PF.band_power(T , λₗ = 1.0 , λᵣ=6.0 ) , 
                PF.∇ₜband_power(T  , bp_value ; λₗ = 1.0 , λᵣ=6.0 ) , 
                PF.∇²ₜband_power(T  , bp_value ; λₗ = 1.0 , λᵣ=6.0 )
                    ) .≈ PF.Dₜband_power(T , λₗ = 1.0 , λᵣ=6.0 , tol =1e-12))


            @test all((PF.band_power(T , λₗ = 1.0 , λᵣ=6.0) , PF.∇ₜband_power(T , λₗ = 1.0 , λᵣ=6.0)) .≈ ChainRulesCore.frule(
                (nt , 1.0), PF.band_power , T ; λₗ = 1.0 , λᵣ=6.0
            ))

            @test all((PF.spectral_ratio(1.0 , 6.0 , T , e_slope = 1.05) , PF.∇ₜspectral_ratio(1.0 , 6.0 , T, e_slope = 1.05)) .≈ ChainRulesCore.frule(
                (nt , nt , nt ,  1.0), PF.spectral_ratio ,1.0 , 6.0 ,  T ; e_slope = 1.05
            ))

        println("ok")

        print("testing weighted average function ...")
                    λ = collect(range(0.8, 2.5, length=500))
                    α = copy(λ)
                    
                    f_test = inv

                    num_integral, _ = quadgk(l -> f_test(l) * PF.ibb(l, T), λ[begin], λ[end], rtol=1e-14)
                    
                    den_integral, _ = quadgk(l -> PF.ibb(l, T), λ[begin], λ[end], rtol=1e-14)
                    den_int_sum = PF.band_power(T , λₗ = λ[begin] , λᵣ=λ[end] , tol =1e-14)
                    @test den_integral ≈ den_int_sum rtol = 1e-8

                    expected = num_integral / den_integral
                                        
                    actual = PF.weighted_average(α, λ, T, PF.ibb , inv)
                    @test actual ≈ expected atol = 1e-5
                   
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
