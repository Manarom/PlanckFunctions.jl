# PlanckFunctions and its derivatives all tempeatures should be in Kelvins, all wavelengths in μm


"""
    PlanckFunctions

A comprehensive suite of tools for evaluating properties of the Planck thermal emission spectrum, 
including spectral intensity (radiance), integrated band power, and their ratios. 

The module provides exact analytical derivatives with respect to both wavelength (λ) 
and temperature (T), optimized for high-performance physics and pyrometry calculations.

### Core Design Principles

* **Derivative Naming Convention**: 
  For any base function `foo`, its first and second temperature derivatives are available via `∇ₜfoo` and `∇²ₜfoo`.
* **Combined Evaluation**: 
  Functions prefixed with `Dₜ` (e.g., `Dₜfoo`) evaluate the base value, the first derivative, and the second derivative simultaneously, returning them as a single `Tuple` to prevent redundant calculations.
* **Functional Operators**: 
  The symbols `∇ₜ` and `∇²ₜ` can be used directly as operators on functions, mapping a base function to its corresponding derivative function (e.g., `PlanckFunctions.∇ₜ(PlanckFunctions.ibb) -> ∇ₜibb`).
  The symbol `∫ₗ` applied on function mapps it to its wavelength integrator , e.g. 
    
    

### Physical Units

* **Wavelength (λ)**: Microns (μm)
* **Temperature (T)**: Kelvins (K)

---

### API Reference

#### Spectral Radiance & Core Functions

* [`ibb`](@ref): Spectral intensity (spectral radiance) in [W/m²⋅sr⋅μm].
* [`∇ₜibb`](@ref): First derivative of spectral intensity with respect to temperature [W/m²⋅sr⋅μm·K].
* [`∇²ₜibb`](@ref): Second derivative of spectral intensity with respect to temperature [W/m²⋅sr⋅μm·K²].

All this functions has also and integrated versions , which can be obtained both by direct name 
or by applying the [`∫ₗ`](@ref) [`∫ₗ`](@ref) operator 

* [`Dₜibb`](@ref): Base functions, its first and second derivative.
* [`ibb!`](@ref),  [`∇ₜibb!`](@ref),  [`∇²ₜibb!`](@ref),  [`Dₜibb!`](@ref) - in-place versions 
* [`∇ₗibb`](@ref): First derivative of spectral intensity with respect to wavelength [W/m²⋅sr⋅μm²].
* [`∇²ₗibb`](@ref): Second derivative of spectral intensity with respect to wavelength [W/m²⋅sr⋅μm³].


#### Integrated Band Power
* [`band_power`](@ref), [`∇ₜband_power`](@ref), [`∇²ₜband_power`](@ref): Total integrated intensity within a specified wavelength region along with its temperature derivatives [W/m²⋅sr].

#### Analytical Ratios (Continuous Pyrometry)
* [`spectral_ratio`](@ref), [`∇ₜspectral_ratio`](@ref), [`∇²ₜspectral_ratio`](@ref): Ratio of spectral intensities at two distinct discrete wavelengths (useful for ratio pyrometers).
* [`spectral_band_ratio`](@ref), [`∇ₜspectral_band_ratio`](@ref), [`∇²ₜspectral_band_ratio`](@ref): Ratio of integrated spectral intensities for two distinct wavelength bands.

#### Averaging & Continuum Mechanics
* [`planck_averaged`](@ref): Planck-averaged spectral quantity (e.g., total emittance).
* [`planck_averaged_attenuation`](@ref): Planck-averaging applied specifically to spectral attenuation coefficients.
* [`rosseland_averaged_attenuation`](@ref): Rosseland-averaging applied to spectral attenuation coefficients.
* [`weighted_average`](@ref): General-purpose function to compute a custom weighted average over a spectrum.

#### Discrete/Tabular Integration (Experimental Data Pyrometry)
* [`planck_weighted`](@ref), [`∇ₜplanck_weighted`](@ref), [`∇²ₜplanck_weighted`](@ref), [`Dₜplanck_weighted`](@ref): Integrates a discrete, wavelength-dependent experimental quantity (such as emissivity or detector response) weighted by the Planck function and its derivatives.
* [`planck_weighted_ratio`](@ref), [`∇ₜplanck_weighted_ratio`](@ref), [`∇²ₜplanck_weighted_ratio`](@ref), [`Dₜplanck_weighted_ratio`](@ref): Computes the ratio of two independent Planck-weighted discrete quantities and its analytical temperature derivatives.
"""
module PlanckFunctions
    using ChainRulesCore
    import Base: /
    export   ∫ibbₗ, Tₖ, 

            ibb , ∇ₜibb, ∇²ₜibb, Dₜibb, # temperature derivatives 

            ibb!, ∇ₜibb!, ∇²ₜibb!, Dₜibb!, # in-place versions 
            
            ∇ₗibb ,  ∇²ₗibb, # wavelength derivatives 

            power , ∇ₜpower ,  ∇²ₜpower , Dₜpower ,  # total intensity and derivatives 

            band_power ,  ∇ₜband_power, ∇²ₜband_power, Dₜband_power, # band integrated Planck function and derivatives 
            
            spectral_ratio, ∇ₜspectral_ratio, ∇²ₜspectral_ratio, Dₜspectral_ratio, # two wavelengths intensities ratio
            
            spectral_band_ratio , ∇ₜspectral_band_ratio, ∇²ₜspectral_band_ratio, Dₜspectral_band_ratio, # two band integral intensities ratio
            
            planck_weighted ,∇ₜplanck_weighted, ∇²ₜplanck_weighted, Dₜplanck_weighted,

            planck_weighted_ratio ,∇ₜplanck_weighted_ratio, ∇²ₜplanck_weighted_ratio,  Dₜplanck_weighted_ratio,

            planck_averaged, planck_averaged_attenuation, rosseland_averaged_attenuation,
            weighted_average,
            weighted_value , ∫ₗ , ∇ₜ , ∇²ₜ  


    const citation = "J.R.Howell,M.P.Menguc,J.R.Howell,M.P.Menguc,K.Daun,R.Siegel. Thermal radiation heat transfer. Seventh edition. 2021 " 
    const citation2 = "Risch, T.K., User's Manual: Routines for Radiative Heat Transfer and Thermometry. NASA/TM-2016-219103. 2016, Edwards, California: Armstrong flight Research Center"
    const citation3 = "https://physics.nist.gov/"
    # FUNDAMENTAL CONSTANTS 
    """
        Planck const
        `h` [$(citation3)] , J Hz^-1
    """
    const h_big  = BigFloat("6.62607015e-34") 
        """
        Speed of light in vacuum
        `c` [$(citation3)] , m s^-1
    """
    const c_big  = BigFloat("299792458.0")
        """
        Boltzmann const 
        `kb` [$(citation3)] , J K^-1    
    """    
    const kB_big = BigFloat("1.380649e-23")    
    const C₁_big   =    2 * h_big * c_big^2 * BigFloat("1e24")
    const C₂_big = (h_big * c_big / kB_big) * BigFloat("1e6")
    const σ_big  = (2 * big(pi)^5 * kB_big^4) / (15 * h_big^3 * c_big^2)

    # Root for base Planck function maximum (ibb): 5 * (1 - exp(-x)) - x = 0
    const x_wien = let
        x = big"4.965114231744276"
        for _ in 1:6
            ex = exp(-x)
            f  = 5 * (one(x) - ex) - x
            df = 5 * ex - one(x)
            x -= f / df
        end
        x
    end

    # Root for 1st temperature derivative peak (∇ₜibb)
    # Equation: 6 * (1 - exp(-x)) - x * (1 + exp(-x)) = 0
    const x_nav1 = let
        x = big"5.421435332832289" # Correct initial guess
        for _ in 1:6
            ex = exp(-x)
            f  = 6 * (one(x) - ex) - x * (one(x) + ex)
            df = 5 * ex - one(x) - x * ex # Analytic first derivative
            x -= f / df
        end
        x
    end
    const x_nav2 = C₂_big /(1000.0 * big"1.957329357329357329357329357329357329357329357329357329357329357329357329357332")

    const C₃_big = C₂_big / x_wien
    const C₄_big = C₁_big / (C₃_big^5 * (exp(x_wien) - one(x_wien)))

    # for (∇ₜibb)
    const C₃_∇₁_big = C₂_big / x_nav1
    const C₄_∇₁_big = (C₁_big * x_nav1 *  exp(x_nav1)) / (C₃_∇₁_big^5 * ( exp(x_nav1) - 1)^2)

    # for (∇²ₜibb)
    const C₃_∇₂_big = C₂_big / x_nav2
    const C₄_∇₂_big = big"1.050264096505020973289772466278833634341988089223353286134188913252405393160527e-10"

    const ħ = Float64(h_big / (2pi)) # J*s
    #const C₁   = 1.191043E8::Float64#(1.191043E8,"W*μm*/m²*sr"," ","Risch","2016"),
    
    """
        C₁ constant for Planck function multiplier in [W⋅μm⁴/(m²⋅sr)]
        source $(citation3) 
    """
    const C₁   = Float64(C₁_big)
    #const C₁   = 2*0.595522001e8 # TRHT v. 1.1.1
    #const C₂ = 14387.752::Float64#(14387.752,"μm**K"," ","Risch","2016"),
    """
        C₂ constant for Planck spectral intensity exponent in [μm⋅K]
        source $(citation3)
    """
    const C₂ = Float64(C₂_big)
    #const C₂ = 14387.7688::Float64# TRHT v. 1.1.1

    """
    C₃ constant of Wien's displacement law [μm⋅K]
    source $(citation3)
""" 
    const C₃ = Float64(C₃_big)
    #const C₃ = 2897.77::Float64#(2897.77,"μm*K"," ","Risch","2016"),
    """
     Constant of the first derivaitve Wien's displacement 
    """
    const C₃_∇₁ = Float64(C₃_∇₁_big)
        """
     Constant of the second derivaitve Wien's displacement 
    """
    const C₃_∇₂ = Float64(C₃_∇₂_big)
    """
        C₄ constant in equation for maximum blackbody intensity [W/(m²⋅μm⋅sr*K⁵)]
        source $(citation3)
    """   
    const C₄ = Float64(C₄_big)
    const C₄_∇₁ = Float64(C₄_∇₁_big)
    const C₄_∇₂ = Float64(C₄_∇₂_big)
    #const C₄ = 4.09567E-12::Float64 #TRHT v. 1.1.1
    #const C₄ = 4.09567E-12::Float64#(4.09567E-12,"W/m^2*μm*sr*K^5"," ", "Risch","2016"),
    """
    Stefan-Boltzmann constant [W/(m²*K⁴)]
source $(citation3) 
"""
    const σ = Float64(σ_big)    
    #const σ  = 5.670367E-8::Float64 v. 1.1.1
    #const σ  = 5.670400E-8::Float64 #(5.670400E-8,"W/(m²*K⁴)"," ", "Risch","2016"), previous constants values 
    const Tₖ = 273.15::Float64 #(273.15,"K"," ", "Risch","2016");
"""
    a₁₂₃(λ::Float64,T::Float64)

# Arguments:

`amat` - matrix of intermediate coefficients size [Nx3]
`λ` - wavelength , μm,  [Nx0]
`T` - temperature,  K

Returns tuple:

(`a1 = C₂/(λ*T)`, `a2 = 1/expm1(a1)`, `a3 = 1 + a2`)
"""
function  a₁₂₃(λ::L , T::F) where {F , L }
    D = promote_type(F , L)
    one_D = one(D)
    a1 = D(C₂/(λ*T))
    a2 = one_D/expm1(a1) # 1/eaxpm1(a)
    a3 = one_D + a2  # exp(a)/expm1(a) = 1 + 1/expm 
    return  (a1 , a2 , a3)
end
    # need branching here ?
    #       T-> 0   :  a1->inf, a2 -> 0 , a3 -> 1 
    #       T -> inf  : a1 -> 0 , a2 -> 1 , a3 -> 1  
    #=if T <= zero(D)
        return (inv(zero(D)), zero(D), one_D)  # T-> 0   :  a1->inf, a2 -> 0 , a3 -> 1 
    end=#
"""
    a₁₂₃!(amat::AbstractMatrix,λ::AbstractVector,T::Number)

In-place filling of the intermediate matrix
`a₁=C₂/(λ*T)`  - amat first column
`a₂ = 1/(eᵃ¹-1)`  - amat second column 
`a₃ = eᵃ¹/(eᵃ¹-1)` - amat third column

# Arguments:

`amat` - matrix of intermediate coefficients size [Nx3],
`λ` - wavelength , μm,  [Nx0],
`T` - temperature,  K,

"""
function a₁₂₃!(amat::AbstractMatrix{D} , λ::AbstractVector , T::Number) where D
    if T <= zero(T)
        @views amat[:, 1] .= inv(zero(D))
        @views amat[:, 2] .= zero(D)
        @views amat[:, 3] .= one(D)
        return amat
    end
    # TRY VIEW WITH broadcating
    a1 = @view amat[: , 1]
    a2 = @view amat[: , 2] 
    a3 = @view amat[: , 3]

    inv_T = inv(T)
    c2_over_T = C₂ * inv_T

    @. begin
        a1 = c2_over_T / λ
        a2 = one(D) / expm1(a1)
        a3 = one(D) + a2
    end
    return amat
end
# BLACKBODY INTENSITY
"""
    ibb(λ,T)

Blackbody spectral intensity (spectral radiance), [W/m²⋅sr⋅μm]
` Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1)` ,
`a₁=C₂/(λ*T)`

# Arguments:

`λ` - wavelength in μm
`T` - temperature in Kelvins
"""
    ibb(λ , T) = ( C₁/expm1( C₂/(λ*T) ) )*λ^-5
    ibb(λ::AbstractVector  , T::Number) = ibb.(λ , T)
    ibb(λ::AbstractVector,T::Base.RefValue{Q}) where Q <: Number = ibb(λ , T[])

"""
    ibb(λ::AbstractVector,amat::AbstractMatrix)

Blackbody spectral intensity (spectral radiance)  with intermediate matrix provided externally, [W/m²⋅sr⋅μm]
`Ibb =  C₁*(λ⁻⁵)*a₂` , 
`a₁=C₂/(λ*T)`  - amat first column, 
`a₂ = 1/(eᵃ¹-1)`  - amat second column 

# Arguments:

`amat` - matrix of intermediate coefficients,  [Nx3]
`λ` - wavelength in μm,  [Nx0]
"""
function ibb(λ::AbstractVector , amat::AbstractMatrix) # internal version with provided coefficients matrix
    a2 = @view amat[: , 2]
    return @. C₁ * a2 * ((1/λ)^5)     
end
        """
        ibb(λ::AbstractVector,T::AbstractVector)

Blackbody spectral intensity (spectral radiance),  [W/m²⋅sr⋅μm]
`Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1)` , where `a₁=C₂/(λ*T)`

# Arguments:

`λ` - wavelength in μm, [Nx0]
`T` - temperature in Kelvins [Mx0]

    """
    ibb(λ::AbstractVector , T::AbstractVector) =  @. (C₁/expm1(C₂/(λ*$transpose(T))))*λ^-5

    
            """
        ibb!(i::AbstractVector , λ::AbstractVector , T::Number)
        
In-place blackbody intensity,  [W/m²⋅sr⋅μm]
`Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1)` , where `a₁=C₂/(λT)`

# Arguments:

`i` - bb intensity vector, [Nx0]
`λ` - wavelength in μm, [Nx0]
`T` - temperature in Kelvins
    """
    ibb!(i::AbstractVector , λ::AbstractVector , T::Number) =  @. i = ibb(λ , T)

    
    """
        ibb!(i::AbstractVector , λ::AbstractVector , amat::AbstractMatrix)
        
In-place blackbody intensity with intermediate coefficients provided externally, [W/m2-sr-mkm]
`Ibb =  C₁*(λ⁻⁵)*a₂` , where
`a₁=C₂/(λ*T)`  - amat first column
`a₂ = 1/(eᵃ¹-1)`  - amat second column 

# Arguments:

`i` - BB intensity, [Nx0]
`λ` - wavelength in μm,  [Nx0]
`amat` - matrix of intermediate coefficients,  [Nx3]
    """
function ibb!(i::AbstractVector , λ::AbstractVector , amat::AbstractMatrix) # this version is used in emissivity approximation 
        a2 = @view amat[: , 2] # Ibb = (λ⁻⁵)* C₁*a₂
        return @. i = C₁ * a2 *((1/λ)^5) 
    end   

    """
            ∇ₜibb(λ,T)

BB intensity first derivative with respect to temperature
dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))
`a₁=C₂/(λ*T)`
`a₂ = 1/(eᵃ¹-1)`   #  1/expm1(a1)
`a₃ = eᵃ¹/(eᵃ¹-1)` #  exp(a)/expm1(a)
dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))

# Arguments:

`λ` - wavelength, μm
`T` - temperature, K
    """
    ∇ₜibb(λ , T) = prod(a₁₂₃(λ , T)) * C₁/(T * λ^5)
    ∇ₜibb(λ::AbstractVector , T::AbstractVector) = @. ∇ₜibb(λ , $transpose(T))    

    """
    ∇ₜibb(λ::AbstractVector,T,amat::AbstractMatrix)

BB intensity first derivative with respect to temperature
with externally provided matrix of intermediate coefficients
`dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))`
`a₁=C₂/(λ*T)`
`a₂ = 1/(eᵃ¹-1)`   
`a₃ = eᵃ¹/(eᵃ¹-1)` 
`dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))`

# Arguments:

`λ` - wavelength in μm, [Nx0]
`T` - temperature in Kelvins
`amat` - matrix of intermediate coefficients, [Nx3]

"""
function ∇ₜibb(λ::AbstractVector , T::Number , amat::AbstractMatrix)
        return C₁*prod(amat; dims=2)./(T*λ.^5)
    end
    """
    ∇ₜibb!(g::AbstractMatrix , λ::AbstractVector , T::AbstractVector)

In-place BB intensity first derivative with respect to temperature
`a₁=C₂/(λ*T)`
`a₂ = 1/(eᵃ¹-1)`   
`a₃ = eᵃ¹/(eᵃ¹-1)` 
`dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))`

# Arguments:

`g` - vector to be filled, [Nx0]
`λ` - wavelength in μm, [Nx0]
`T` - temperature in Kelvins     
"""
function ∇ₜibb!(g::AbstractMatrix , λ::AbstractVector , T::AbstractVector)
        # instance version of Planck function first derivative with respect to T
        size(g) == (length(λ), length(T)) || throw(DimensionMismatch())
        @inbounds for j in eachindex(T)
            for i  in eachindex(λ) 
                g[i , j] = ∇ₜibb(λ[i] , T[j])
            end
        end   
    end
    ∇ₜibb!(g::AbstractVector , λ::AbstractVector , T) = @. g = ∇ₜibb(λ , T)
    # this version uses Vector for T because in this way the handle to the optimization varibale is implemented
    # T should be one-element array!
    """
    ∇ₜibb!(g::AbstractVector , λ::AbstractVector , T , amat::AbstractMatrix)

In-place bb intensity first derivative with respect to temperature
with externally provided amat  - matrix with columns a₁,a₂,a₃

`dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))`,
`a₁=C₂/(λ*T)` ,
`a₂ = 1/(eᵃ¹-1)`   ,
`a₃ = eᵃ¹/(eᵃ¹-1)` ,
`dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))`,

# Arguments:
g - to be filled, [Nx0]
λ - wavelength in μm, [Nx0]
T - temperature in Kelvins
amat - matrix of intermediate coefficients, [Nx3]
"""
function ∇ₜibb!(g::AbstractVector , λ::AbstractVector , T , amat::AbstractMatrix)# for fixed value of temperature
        # instance version of Planck function first derivative with respect to T
        #a = zeros(3)
        #t = T[1]
        prod!(g , amat) # this puts all amat rows with column-wise product in g
        @. g *= C₁/(T * λ^5) # a₃*a₁*a₂*C₁*(1/λ⁵)*(1/T)
        return g
    end
    """
    ∇ₜibb!(g::AbstractVector , T , amat::AbstractMatrix , bb_intensity::AbstractVector)

In-place bb intensity first derivative with respect to temperature
with externally provided both `amat`  - matrix with columns a₁,a₂,a₃
and `bb_intensity`

dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))
a₁=C₂/(λ*T)
a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))
as far as Ibb = C₁*a₂/λ⁵
dIbb/dT = a₃*a₁*C₁*(a₂/λ⁵)*(1/T)=a₃*a₁*Ibb/T

# Arguments:
g - to be filled, [Nx0]
λ - wavelength in μm, [Nx0]
T - temperature in Kelvins
amat - matrix of intermediate coefficients, [Nx3]
"""
function ∇ₜibb!(g::AbstractVector , T , amat::AbstractMatrix , bb_intensity::AbstractVector)
        #t = T[1]
        a1 = view(amat,:,1)
        a3 = view(amat,:,3)
        @. g = a1 * a3 * bb_intensity/T  # dIbb/dT = a₃*a₁*Ibb/T
        return g
    end
"""
    ∇²ₜibb(λ,T)

BB intensity second derivative with respect to temperature

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]` ,

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]`,
`a₁=C₂/(λ*T)`,
`a₂ = 1/(eᵃ¹-1)`,  
`a₃ = eᵃ¹/(eᵃ¹-1)` ,

`d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]`

# Arguments :
λ - wavelength in μm
T - tmperature in Kelvins         
"""
    function ∇²ₜibb(λ,T)
        (a1 , a2 , a3) = a₁₂₃(λ,T)#ibb = C₁*a[2]*(l^-5)
        return (a1 * (2.0 * a3 - 1.0) - 2.0) * a1 * a2 * a3 * C₁ / ((T^2) * λ^5)
    end
    ∇²ₜibb(λ::AbstractVector , T::AbstractVector) = @. ∇²ₜibb(λ , $transpose(T)) 
    """
    ∇²ₜibb!(h::AbstractVector , λ::AbstractVector , T::Number)

In-place bb intensity second order derivative with respect to temperature

# Arguments :
h  - to be filled, [Nx0]
λ - wavelength in μm, [Nx0]
T - tmperature in Kelvins  
"""
∇²ₜibb!(h::AbstractVector , λ::AbstractVector , T::Number)  = @. h =  ∇²ₜibb(λ,T)# secpnd derivative for the fixed value of temperature

    """
    ∇²ₜibb!(h::AbstractMatrix{Float64} , λ::AbstractVector{Float64}, T::AbstractVector{Float64})

In-place bb intensity second order derivative with respect to temperature

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]`,

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]`,
`a₁=C₂/(λ*T)`,
`a₂ = 1/(eᵃ¹-1) `  ,
`a₃ = eᵃ¹/(eᵃ¹-1)` ,

`d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]`

# Arguments :
`h`  - to be filled, [Nx0]
`λ` - wavelength in μm, [Nx0]
`T `- tmperature in Kelvins   
"""
function ∇²ₜibb!(h::AbstractMatrix , λ::AbstractVector , T::AbstractVector)
        size(h) == (length(λ), length(T)) || throw(DimensionMismatch())
        @inbounds for j in eachindex(T)
            for i  in eachindex(λ)         
                h[i,j] = ∇²ₜibb(λ[i] , T[j]) 
            end
        end 
        return h    
    end
    """
    ∇²ₜibb!(h::AbstractVector{Float64}, λ::AbstractVector{Float64} , T::Float64 ,amat::AbstractMatrix{Float64})::Nothing

In-place bb intensity second order derivative with respect to temperature with 
intermediate matrix provided externally

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]`,

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]`,
`a₁=C₂/(λ*T)`,
`a₂ = 1/(eᵃ¹-1)`   ,
`a₃ = eᵃ¹/(eᵃ¹-1)` ,

`d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]`

# Arguments :

`h`  - to be filled, [Nx0]
`λ `- wavelength in μm, [Nx0]
`T` - temperature in Kelvins
`amat` - matrix of intermediate coefficients,  [Nx3]
"""
function ∇²ₜibb!(h::AbstractVector , λ::AbstractVector , T::Number , amat::AbstractMatrix)
            (a1, a2, a3) = ntuple(i -> view(amat, :, i), 3)
            @. h = C₁ * (a1 * a2 * a3) *  ( a1 * (2a3 - 1.0) - 2.0 ) /( (T^2)*λ.^5 ) # C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]
        return h        
    end
    """
    ∇²ₜibb!(h::AbstractVector{Float64},T::Float64,amat::AbstractMatrix{Float64},∇i::AbstractVector{Float64})::Nothing

In-place bb intensity second order derivative with respect to temperature 
with provided both the intermediate matrix amat and the the Planck function first derivative

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]`

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]`
`a₁=C₂/(λ*T)`
`a₂ = 1/(eᵃ¹-1)`   ,
`a₃ = eᵃ¹/(eᵃ¹-1)` ,
`d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]`
as far as 
    `Ibb = (λ⁻⁵)* C₁*a₂`
and 
    `dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T)) = a₃*a₁*Ibb/T` 
hense

`    d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
        = [a₃*a₁*Ibb/T^2]*[a₁*(2*a₃ - 1))-2] 
            = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] `

# Arguments :

`h`  - to be filled, [Nx0]
`λ` - wavelength in μm, [Nx0]
`amat` - matrix of intermediate coefficients,  [Nx3]
`∇i` - vector of bb intensity first derivatives, [Nx0]
"""
function ∇²ₜibb!(        h::AbstractVector , 
                        T::Number , 
                        amat::AbstractMatrix,
                        ∇i::AbstractVector)
        # instance version of Planck function second derivative with respect to T
        # with supplied coefficients matrix
            a1 = view(amat,:,1)
            a3 = view(amat,:,3)
            @. h = (a1 * (2.0*a3 -1.0) -2.0) * ∇i/T  # h = (∇Ibb/T)*[a₁*(2*a₃ - 1))-2] 
        return h       
    end
    """
    rosseland_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number)

Evaluates the Rosseland-averaged spectral attenuation coefficient (the summation of
spectral scattering and absorption coefficients) `α(λ)` for temperature `T`:

`αᵣ = (∫(1/α(λ))⋅∇ₜibb(λ,T)dλ/∫∇ₜibb(λ,T)dλ)⁻¹`  
"""
rosseland_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number) = inv(weighted_average(α , λ , T , ∇ₜibb , inv))

"""
    planck_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number)

Planck-averaged spectral attenuation coefficient (the summation of
spectral scattering and absorption coefficients) α(λ) for temperature T:

`αᵣ = (∫(1/α(λ))⋅ibb(λ,T)dλ/∫ibb(λ,T)dλ)⁻¹`   
"""
planck_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number) = inv(weighted_average(α , λ , T , ibb , inv)) #double inversion

"""
    planck_averaged(x::AbstractVector, λ::AbstractVector,T::Number)

Evaluates the Planck-averaged value of `x(λ)` for temperature `T`:

`xᵣ = ∫x(λ)ibb(λ,T)dλ/∫ibb(λ,T)dλ `

E.g. can be used to evaluate the integral  from  spectral emissivity
"""
planck_averaged(x::AbstractVector, λ::AbstractVector,T::Number) = weighted_average(x , λ , T , ibb , identity)
const Dfunctions = Union{typeof(ibb),typeof(∇ₜibb),typeof(∇²ₜibb)}
"""
    weighted_average(α::AbstractVector, 
                        λ::AbstractVector,
                        T::Number, 
                        g::Dfunctions,
                        f::F = identity) where F

Generic function to evaluate the averaged value of some `f(x)` function of variable `x` dependent
on wavelength `λ` for temperature `T`. 
Uses linear approximation for the discrete variable and square polynomial 
for the g function 

`xᵣ = ∫f(x)g(λ,T)dλ/∫g(λ,T)dλ`
the default value of f is `identity`, e.g. if f = inv:
`xᵣ = ∫g(λ,T)/x(λ)dλ/∫g(λ,T)dλ`

"""
function weighted_average(α , λ , T , g , f = identity)
        (s , sn) = weighted_value(α , λ ,T , g , f)
        return s/sn
end

"""
    weighted_value(     α::AbstractVector, 
                        λ::AbstractVector,
                        T , 
                        g::Union{typeof(ibb),typeof(∇ₜibb),typeof(∇²ₜibb)},
                        f::F = identity) where F <: Function

Returns the tuple of  `(∫f(α)g(λ,T)dλ , ∫g(λ,T)dλ)`

"""
function weighted_value(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D, 
                        g::Dfunctions,
                        f::F = identity) where {F   , D <: Number , Q <: Number , L <: Number}

    N = length(α)
    @assert length(λ) == N

    DD = promote_type(Q , D , L)

    s = zero(DD)
    sn = zero(DD)

    (lmin , lmax) = extrema(λ)
    (_ , norm_val) = maximum_on_interval(g , lmin , lmax , T)    
    nrm = inv(norm_val)

    g_start = g(λ[begin], T) * nrm
    
    inv2 = inv(DD(2))
    inv3 = inv(DD(3))
    inv4 = inv(DD(4))

    @inbounds for i in 1:N-1
        lstart = λ[i]
        lend = λ[i + 1]
        
        h = lend - lstart
        h_half = h * inv2

        b1 = DD(f(α[i]))
        b2 = (DD(f(α[i+1])) - b1) / h

        g_end = g(lend, T) * nrm
        g_centre = g(lstart + h_half, T) * nrm
        
        (a1, a2, a3) = second_order_polynomial_fit(
            zero(D), h_half, h,
            g_start, g_centre, g_end
        )

        g_start = g_end                                        

        c2 = a1 * b1
        c3 = (a1 * b2 + a2 * b1) * inv2
        c4 = (a2 * b2 + a3 * b1) * inv3
        c5 = a3 * b2 * inv4

        s  += @evalpoly(h, zero(DD), c2, c3, c4, c5)
        sn += @evalpoly(h, zero(DD), a1, a2 * inv2, a3 * inv3)
    end 
    return (s * norm_val , sn * norm_val)
end
"""
    weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D, 
                        g::Dfunctions) where { D <: Number , Q <: Number , L <: Number}

  `∫α⋅g(λ,T)dλ` where g ∈ {ibb , ∇ₜibb , ∇²ₜibb}
"""
function weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D, 
                        g::Dfunctions) where { D <: Number , Q <: Number , L <: Number}

    N = length(α)
    @assert length(λ) == N

    DD = promote_type(Q , D , L)
    s = zero(DD)
    (lmin , lmax) = extrema(λ)
    (_ , norm_val) = maximum_on_interval(g , lmin , lmax , T)    
    nrm = inv(norm_val)

    g_start = g(λ[begin], T) * nrm
    
    inv2 = inv(DD(2))
    inv3 = inv(DD(3))
    inv4 = inv(DD(4))

    @inbounds for i in 1:N-1
        lstart = λ[i]
        lend = λ[i + 1]
        
        h = lend - lstart
        h_half = h * inv2

        b1 = DD(α[i])
        b2 = (DD(α[i+1]) - b1) / h

        g_end = g(lend, T) * nrm
        g_centre = g(lstart + h_half, T) * nrm
        
        (a1, a2, a3) = second_order_polynomial_fit(
            zero(D), h_half, h,
            g_start, g_centre, g_end
        )

        g_start = g_end                                        
        c2 = a1 * b1
        c3 = (a1 * b2 + a2 * b1) * inv2
        c4 = (a2 * b2 + a3 * b1) * inv3
        c5 = a3 * b2 * inv4
        s  += @evalpoly(h, zero(DD), c2, c3, c4, c5)
    end 
    return s * norm_val 
end



"""
    planck_weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}
Evaluates `∫α(λ)⋅ibb(λ,T)dλ` using numeric integration

"""
planck_weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}  = weighted(α , λ , T , ibb)
 """
    planck_weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}\

Evaluates `∫α(λ)⋅∇ₜibb(λ,T)dλ` using numeric integration
"""
∇ₜplanck_weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}  = weighted(α , λ , T , ∇ₜibb)           
"""
    ∇²ₜplanck_weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}

Evaluates `∫α(λ)⋅∇²ₜibb(λ,T)dλ` using numeric integration
"""
∇²ₜplanck_weighted(α::AbstractVector{Q}, 
                        λ::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}  = weighted(α , λ , T , ∇²ₜibb)                                      
"""
    weighted_values_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D, 
                        g::Dfunctions,
                        f::F = identity) where {F   , D <: Number , Q <: Number , L <: Number}


Returns the ratio of two weighted values  `(∫f(α1)g(λ,T)dλ1 / ∫f(α2)g(λ,T)dλ2)`
"""
function weighted_values_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D, 
                        g::Dfunctions,
                        f::F = identity) where {F   , D <: Number , Q <: Number , L <: Number}

    return first(weighted_value(α1 , λ1 , T , g , f))/first(weighted_value(α2 , λ2 , T , g , f))

end
"""
    planck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}

Returns the ratio of two weighted values  `(∫α1(λ)ibb(λ,T)dλ1 / ∫α2(λ)ibb(λ,T)dλ2)`
"""
function planck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}
    return planck_weighted(α1 , λ1 , T )/planck_weighted(α2 , λ2 , T)
end

"""
    ∇ₜplanck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}

Returns the ratio of two weighted values  `d/dT [ (∫α1(λ)ibb(λ,T)dλ1 / ∫α2(λ)ibb(λ,T)dλ2) ]`
"""
function ∇ₜplanck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}

    (f1 , df1 , _ ) = Dₜplanck_weighted(α1 , λ1 , T)
    (f2 , df2 , _ ) = Dₜplanck_weighted(α2 , λ2 , T)

    return _spectral_ratio_first_derivative(f1 , df1 , f2 , df2)
end
"""
    ∇²ₜplanck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}

Returns the ratio of two weighted values  `d²/dT² [ (∫α1(λ)ibb(λ,T)dλ1 / ∫α2(λ)ibb(λ,T)dλ2) ]`
"""
function ∇²ₜplanck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D) where { D <: Number , Q <: Number , L <: Number}

    (f1 , df1 , d2f1 ) = Dₜplanck_weighted(α1 , λ1 , T)
    (f2 , df2 , d2f2 ) = Dₜplanck_weighted(α2 , λ2 , T)
    return _spectral_ratio_second_derivative(f1 , df1 , d2f1 , f2 , df2 , d2f2)
end
"""
    Dₜplanck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D ) where { D <: Number , Q <: Number , L <: Number}

Returns the tuple of the value, first and second derivatives of two weighted values ratio :

# Arguments 
`α1` - first wavelength dependent quantity, 
`λ1` - first wavelengths range, μm
`α2`-  second wavelength dependent quantity, 
`λ2`-  second wavelengths range, μm
`T` - temperature , K

# Returns
``( ∫α1⋅g(λ,T)dλ1 / ∫f(α2)g(λ,T)dλ2 , 

   d/dT [∫α1⋅g(λ,T)dλ1/∫α2⋅g(λ,T)dλ2] , 

   d²/dT² [ ∫α1⋅g(λ,T)dλ1/∫α2⋅g(λ,T)dλ2 ]
)``
"""
function Dₜplanck_weighted_ratio(
                        α1::AbstractVector{Q}, 
                        λ1::AbstractVector{L},
                        α2::AbstractVector{Q}, 
                        λ2::AbstractVector{L},
                        T::D ) where { D <: Number , Q <: Number , L <: Number}

    (f1 , df1 , d2f1 ) = Dₜplanck_weighted(α1 , λ1 , T)
    (f2 , df2 , d2f2 ) = Dₜplanck_weighted(α2 , λ2 , T)

    return    return (    
                f1/f2 ,
                _spectral_ratio_first_derivative(f1 , df1 , f2 , df2) , 
                _spectral_ratio_second_derivative(f1, df1, d2f1,
                                                            f2, df2, d2f2)
        )

end
"""
    Dₜplanck_weighted_values(e::AbstractVector{Q}, λ::AbstractVector{L}, T::D) where {Q <: Number, L <: Number, D <: Number}

Evaluates the integral of `e(λ)` over three weighting functions of 
Planck's function ,  its first- and second- derivaitve: 
`ibb`, `∇ₜibb` и `∇²ₜibb`

Returns six-elements tuple (integrals and normalizers) :
`( ∫α(λ)ibb(λ,T)dλ , ∫α(λ)∇ₜibb(λ,T)dλ , ∫α(λ)∇²ₜibb(λ,T)dλ , 
   ∫ibb(λ,T)dλ , ∫∇ₜibb(λ,T)dλ , ∫∇²ₜibb(λ,T)dλ)`
"""
function Dₜplanck_weighted_values(e::AbstractVector{Q}, λ::AbstractVector{L}, T::D) where {Q <: Number, L <: Number, D <: Number}
    N = length(e)
    @assert length(λ) == N

    DD = promote_type(Q, D, L)

    # starting summations
    s1, sn1 = zero(DD), zero(DD) #  ibb
    s2, sn2 = zero(DD), zero(DD) #  ∇ₜibb
    s3, sn3 = zero(DD), zero(DD) #  ∇²ₜibb

    # normalizators 
    (lmin, lmax) = extrema(λ)
    (_, norm_val1) = maximum_on_interval(ibb, lmin, lmax, T)    
    (_, norm_val2) = maximum_on_interval(∇ₜibb, lmin, lmax, T)    
    (_, norm_val3) = maximum_on_interval(∇²ₜibb, lmin, lmax, T)    

    nrm1 = inv(norm_val1)
    nrm2 = inv(norm_val2)
    nrm3 = inv(norm_val3)

    (g1_start , g2_start  , g3_start) = Dₜibb(λ[begin], T) 
    g1_start *= nrm1
    g2_start *= nrm2
    g3_start *= nrm3

    inv2 = inv(DD(2))
    inv3 = inv(DD(3))
    inv4 = inv(DD(4))

    @inbounds for i in 1:N-1
        lstart = λ[i]
        lend = λ[i + 1]
        
        h = lend - lstart
        h_half = h * inv2

        # e approximaiton coefficients 
        b1 = DD(e[i])
        b2 = (DD(e[i+1]) - b1) / h

        # functions values 
        (g1_end , g2_end , g3_end) = Dₜibb(lend , T)
        g1_end *= nrm1
        g2_end *= nrm2
        g3_end *= nrm3

        g1_centre , g2_centre , g3_centre =  Dₜibb(lstart + h_half , T)
        g1_centre *= nrm1
        g2_centre *= nrm2
        g3_centre *= nrm3
        
        # hardcoded polyfit
        (a1_1, a2_1, a3_1) = second_order_polynomial_fit(zero(D), h_half, h, g1_start, g1_centre, g1_end)
        (a1_2, a2_2, a3_2) = second_order_polynomial_fit(zero(D), h_half, h, g2_start, g2_centre, g2_end)
        (a1_3, a2_3, a3_3) = second_order_polynomial_fit(zero(D), h_half, h, g3_start, g3_centre, g3_end)

        # returns back 
        g1_start = g1_end                                        
        g2_start = g2_end                                        
        g3_start = g3_end                                        

        # ibb
        c2_1 = a1_1 * b1
        c3_1 = (a1_1 * b2 + a2_1 * b1) * inv2
        c4_1 = (a2_1 * b2 + a3_1 * b1) * inv3
        c5_1 = a3_1 * b2 * inv4

        # ∇ₜibb
        c2_2 = a1_2 * b1
        c3_2 = (a1_2 * b2 + a2_2 * b1) * inv2
        c4_2 = (a2_2 * b2 + a3_2 * b1) * inv3
        c5_2 = a3_2 * b2 * inv4

        # ∇²ₜibb
        c2_3 = a1_3 * b1
        c3_3 = (a1_3 * b2 + a2_3 * b1) * inv2
        c4_3 = (a2_3 * b2 + a3_3 * b1) * inv3
        c5_3 = a3_3 * b2 * inv4

        s1  += @evalpoly(h, zero(DD), c2_1, c3_1, c4_1, c5_1)
        sn1 += @evalpoly(h, zero(DD), a1_1, a2_1 * inv2, a3_1 * inv3)
        s2  += @evalpoly(h, zero(DD), c2_2, c3_2, c4_2, c5_2)
        sn2 += @evalpoly(h, zero(DD), a1_2, a2_2 * inv2, a3_2 * inv3)
        s3  += @evalpoly(h, zero(DD), c2_3, c3_3, c4_3, c5_3)
        sn3 += @evalpoly(h, zero(DD), a1_3, a2_3 * inv2, a3_3 * inv3)
    end 
    
    return (
        s1 * norm_val1, s2 * norm_val2, s3 * norm_val3, 
        sn1 * norm_val1, sn2 * norm_val2, sn3 * norm_val3
    )
end

"""
    Dₜplanck_weighted(e::AbstractVector{Q}, λ::AbstractVector{L}, T::D) where {Q <: Number, L <: Number, D <: Number}

Evaluates the integral of `e(λ)` over three weighting functions of 
Planck's function ,  its first- and second- derivaitve: 
`ibb`, `∇ₜibb` и `∇²ₜibb`

Returns six-elements tuple (integrals and normalizers) :
`( ∫α(λ)ibb(λ,T)dλ , ∫α(λ)∇ₜibb(λ,T)dλ , ∫α(λ)∇²ₜibb(λ,T)dλ)`
"""
function Dₜplanck_weighted(e::AbstractVector{Q}, λ::AbstractVector{L}, T::D) where {Q <: Number, L <: Number, D <: Number}
    N = length(e)
    @assert length(λ) == N

    DD = promote_type(Q, D, L)

    # starting summations
    s1 = zero(DD) #  ibb
    s2 = zero(DD) #  ∇ₜibb
    s3 = zero(DD) #  ∇²ₜibb

    # normalizators 
    (lmin, lmax) = extrema(λ)
    (_, norm_val1) = maximum_on_interval(ibb, lmin, lmax, T)    
    (_, norm_val2) = maximum_on_interval(∇ₜibb, lmin, lmax, T)    
    (_, norm_val3) = maximum_on_interval(∇²ₜibb, lmin, lmax, T)    

    nrm1 = inv(norm_val1)
    nrm2 = inv(norm_val2)
    nrm3 = inv(norm_val3)

    (g1_start , g2_start  , g3_start) = Dₜibb(λ[begin], T) 

    g1_start *= nrm1
    g2_start *= nrm2
    g3_start *= nrm3

    inv2 = inv(DD(2))
    inv3 = inv(DD(3))
    inv4 = inv(DD(4))

    @inbounds for i in 1:N-1
        lstart = λ[i]
        lend = λ[i + 1]
        
        h = lend - lstart
        h_half = h * inv2

        # e approximaiton coefficients 
        b1 = DD(e[i])
        b2 = (DD(e[i+1]) - b1) / h

        # functions values 
        (g1_end , g2_end , g3_end) = Dₜibb(lend , T)

        g1_end *= nrm1
        g2_end *= nrm2
        g3_end *= nrm3

        g1_centre , g2_centre , g3_centre =  Dₜibb(lstart + h_half , T)
        g1_centre *= nrm1
        g2_centre *= nrm2
        g3_centre *= nrm3
        
        # hardcoded polyfit
        (a1_1, a2_1, a3_1) = second_order_polynomial_fit(zero(D), h_half, h, g1_start, g1_centre, g1_end)
        (a1_2, a2_2, a3_2) = second_order_polynomial_fit(zero(D), h_half, h, g2_start, g2_centre, g2_end)
        (a1_3, a2_3, a3_3) = second_order_polynomial_fit(zero(D), h_half, h, g3_start, g3_centre, g3_end)

        # returns back 
        g1_start = g1_end                                        
        g2_start = g2_end                                        
        g3_start = g3_end                                        

        # ibb
        c2_1 = a1_1 * b1
        c3_1 = (a1_1 * b2 + a2_1 * b1) * inv2
        c4_1 = (a2_1 * b2 + a3_1 * b1) * inv3
        c5_1 = a3_1 * b2 * inv4

        # ∇ₜibb
        c2_2 = a1_2 * b1
        c3_2 = (a1_2 * b2 + a2_2 * b1) * inv2
        c4_2 = (a2_2 * b2 + a3_2 * b1) * inv3
        c5_2 = a3_2 * b2 * inv4

        # ∇²ₜibb
        c2_3 = a1_3 * b1
        c3_3 = (a1_3 * b2 + a2_3 * b1) * inv2
        c4_3 = (a2_3 * b2 + a3_3 * b1) * inv3
        c5_3 = a3_3 * b2 * inv4

        s1  += @evalpoly(h, zero(DD), c2_1, c3_1, c4_1, c5_1)
        s2  += @evalpoly(h, zero(DD), c2_2, c3_2, c4_2, c5_2)
        s3  += @evalpoly(h, zero(DD), c2_3, c3_3, c4_3, c5_3)
    end 
    
    return (
        s1 * norm_val1, s2 * norm_val2, s3 * norm_val3
    )
end


function maximum_on_interval(g::Dfunctions , lmin , lmax , T)
    λmax = λₘ(g , T) # selected function maximum value   
    
    Lmax = if λmax <= lmax
                if λmax >= lmin
                    λmax
                else # λmax < lmin
                    lmin    
                end
            else
                lmax
            end
    return (Lmax , g(Lmax , T))
end
    """
    λₘ(T::Number)

The wavelength (in μm) of bb intensity maximum vs temperature `T` 
`argmax(Planck(T))`  - Wien's displacement law

# Arguments:

`T` - temperature in Kelvins
"""
λₘ(T::Number) = C₃/T
"""
    λₘ(::Union{typeof(ibb) , typeof(∇ₜibb) , typeof(∇²ₜibb)} , T::Number)

The wavelength (in μm) of bb intensity  , its first or second derivaitve maximum vs temperature `T` 

# Arguments:

`T` - temperature in Kelvins
"""
function λₘ(::Dfunctions , T::Number)  end

λₘ(::typeof(ibb) , T::Number) = C₃/T
λₘ(::typeof(∇ₜibb) , T::Number) = C₃_∇₁/T
λₘ(::typeof(∇²ₜibb) , T::Number) = C₃_∇₂/T

    """
    bb_max(T::Number)

Blackbody instensity at maximum value 
"""
bb_max(T::Number) = C₄ * T^5
"""
    bb_max(bb_function , T)

Returns the value of function at a maximum 
`bb_function` can be `ibb` ,   `∇ₜibb` ,  `∇²ₜibb`

"""
bb_max(::typeof(ibb) , T::Number) =  C₄ * T^5
bb_max(::typeof(∇ₜibb), T::Number) = C₄_∇₁ * T^4
bb_max(::typeof(∇²ₜibb) , T::Number) = C₄_∇₂ * T^3

    """
    tₘ(λ)

The temperature of BB having maximum at wavelength `λ` in Kelvins
"""
    function tₘ(λ)
        # temperature (in Kelvins) of BB with intencity maximum at λ μm  
        C₃./λ
    end
        """
        temperature(i)

The bb temperature for intensity `i` corresponding to 
wavelength `λ` 
    """
    function temperature(i , λ)
        return C₂ / (λ * log1p(C₁ / (i * λ^5)))
    end
    """
    bright_temperature(i, λ; ϵ=1.0)

Evaluates bright temperature for single wavelength pyrometer 
`i` - measured intensity ,
`λ` - wavelength, μm 
`ϵ` - emissivity 
    
"""
function bright_temperature(i , λ ; ϵ::Number = 1.0)
        return C₂ / (λ * log1p(ϵ * C₁ / (i * λ^5)))
    end
"""
    ∇ₗibb(λ,T)


BB intensity first derivative with respect to the wavelength

# Arguments:
`λ` - wavelength, μm
`T` - temperature, K
"""
    function ∇ₗibb(λ,T)
        (a1 , a2 , a3) = a₁₂₃(λ,T)
        i = C₁ * a2 * (λ^-5)
        return (i/λ) * (a1 * a3 - 5.0)#
    end
    """
    ∇²ₗibb(λ,T)

BB intensity second derivative with respect to the wavelength

# Arguments:

`λ` - wavelength, μm
`T` - temperature, K
"""
function  ∇²ₗibb(λ,T)
        # second derivative of Planck function with respect to wavelength
        #local a,e2,e3
        (a1 , a2 , a3) = a₁₂₃(λ,T)
        return C₁ * a2 * (a1 * a3*(2a1 * a3 - a1 - 12.0) + 30.0)/(λ^7)
    end
    """
    Dₗibb(λ,T)

Returns a three-element tuple of (1.bb intensity,2.its first and 3.second derivative 
with respect to the wavelentgh)

# Arguments:
`λ` - wavelength, μm
`T` - temperature, K
"""
function Dₗibb(λ , T)
        (a1 , a2 , a3) = a₁₂₃(λ,T)
        i = C₁ * a2 *((1/λ)^5)
        return (
            i,  # Planck function
            (i/λ) * (a1 * a3 - 5.0), # first derivative
            C₁ * a2 * (a1 * a3*(2a1 * a3 - a1 - 12.0) + 30.0)/(λ^7) # second derivative
        )
    end
    function Dₗibb(λ::AbstractVector,T::AbstractVector{Q}) where Q
        # returns spectral intensity and its first and second derivatives with respect to the wavelength
        N , M = length(λ) , length(T)
        i , d1i , d2i = ntuple( _ -> Matrix{Q}(undef , N , M) , 3)

        @inbounds for jjj in eachindex(T)
            for iii in eachindex(λ) 
                l = λ[iii]
                t = T[jjj]
                (i[iii,jjj] , d1i[iii,jjj] , d2i[iii,jjj]) =  Dₗibb(l , t)
            end
        end

        return (i , d1i , d2i)
    end
    """
    power(T)

Integral (over the wavelength) intensity of BB (radiance)  at temperature T

Units: W/(m²⋅sr)

# Arguments:
`T` - temperature, K
"""
power(T) = σ*(T^4)/π

    """
    ∇ₜpower(T)

Total intensity first
derivative of BB (radiance)  at temperature T

Units: W/(m²⋅sr⋅K)

# Arguments:

`T` - temperature, K

"""
    ∇ₜpower(T) = 4*σ*(T^3)/π

    """
    ∇²ₜpower(T)

Total  intensity second 
derivative of BB (radiance)  at temperature T

Units: W/(m²⋅sr⋅K)

# Arguments:
`T` - temperature, K
"""
∇²ₜpower(T) = 12.0 * σ * (T^2) / pi

"""
    Dₜpower(T)

Tuple of total intensity , fisrt and second derivative with respect to temeperature
"""
Dₜpower(T) = (power(T) , ∇ₜpower(T)   , ∇²ₜpower(T) )
    # dummy type to drop setindex
    struct Drop end
    @inline Base.setindex!(::Drop, _ , idx...) = nothing

    function _Dₜibb_core!(out1, out2, out3, λ::AbstractVector, T::Number)
        N = length(λ)
        # Проверяем размеры только у реальных векторов
        !(out1 isa Drop) && (length(out1) != N) && throw(DimensionMismatch())
        !(out2 isa Drop) && (length(out2) != N) && throw(DimensionMismatch())
        !(out3 isa Drop) && (length(out3) != N) && throw(DimensionMismatch())

        @inbounds for iii in eachindex(λ) 
            l = λ[iii]
            a1, a2, a3 = a₁₂₃(l, T)
            _i = C₁ * a2 * (l^-5)
            out1[iii] = _i
            _i2 = a1 * a3 * _i / T
            out2[iii] = _i2
            out3[iii] = (a1 * (2.0 * a3 - 1.0) - 2.0) * _i2 / T
        end
        return (out1 , out2 , out3)
    end

   """
    Dₜibb( λ::Number, T::Number)

returns a tuple of (value, first derivative, second derivative)
"""
function Dₜibb( λ::Number, T::Number)
            a1, a2, a3 = a₁₂₃(λ, T)
            _i = C₁ * a2 * (λ ^-5)
            out1 = _i
            _i2 = a1 * a3 * _i / T
            out2 = _i2
            out3 = (a1 * (2.0 * a3 - 1.0) - 2.0) * _i2 / T
        return (out1 , out2 , out3)
    end
       """
    Dₜibb( λ::Number, T::Number , skip_second_derivative::Val{true})

Skips the second derivative evaluation `Dₜibb( λ, T , Val(true))`
"""
function Dₜibb( λ::Number, T::Number , skip_second_derivative::Val{true})
            a1, a2, a3 = a₁₂₃(λ, T)
            _i = C₁ * a2 * (λ ^-5)
            out1 = _i
            _i2 = a1 * a3 * _i / T
            out2 = _i2
        return (out1 , out2 , nothing)
    end
        """
        Dₜibb!(input_tuple, λ::AbstractVector,T)

In-place filling the tuple of (bb intensity, its first ,and second ) derivatives with 
respect to temperature

# Arguments:
`input_tuple`, {Nx0 vector or nothing , Nx0 vector or nothing, Nx0 vector or nothing}
`λ` - wavelength, μm, [Nx0]
`T` - temperature, K   
Out:
`input_tuple` filled with ibb , dIbb/dT , d²Ibb/dT²

This version is slightly faster than calling each derivative separately
    """ 
    Dₜibb!(input_tuple::Tuple{AbstractVector , AbstractVector , AbstractVector}, 
                    λ::AbstractVector , 
                    T::Number) = _Dₜibb_core!(input_tuple[1],input_tuple[2],input_tuple[3] , λ , T)

    Dₜibb!(input_tuple:: Tuple{Nothing , AbstractVector , AbstractVector}, 
                    λ::AbstractVector , T) = _Dₜibb_core!(Drop() , input_tuple[2] ,input_tuple[3] , λ , T)

    Dₜibb!(input_tuple:: Tuple{AbstractVector , Nothing , AbstractVector}, 
        λ::AbstractVector , T) =_Dₜibb_core!(input_tuple[1] , Drop() ,input_tuple[3] , λ , T)

    Dₜibb!(input_tuple:: Tuple{ AbstractVector , AbstractVector , Nothing}, λ::AbstractVector , T) = _Dₜibb_core!( input_tuple[1] ,input_tuple[2] , Drop() , λ , T)

    function Dₜibb!(    input_tuple::Tuple{AbstractMatrix,AbstractMatrix,AbstractMatrix} , 
                        λ::AbstractVector , 
                        T::AbstractVector)
        N = length(λ)
        M = length(T)
        for m in input_tuple 
            (size(m) != (N , M)) && throw(DimensionMismatch())
        end
        @inbounds for jjj in eachindex(T)
            t = T[jjj]
            i1  = @view input_tuple[1][:,jjj]
            i2  = @view input_tuple[2][:,jjj]
            i3  = @view input_tuple[3][:,jjj]
            _Dₜibb_core!(i1 , i2 , i3 , λ , t)
        end
        return input_tuple
    end
    function Dₜibb(λ::AbstractVector , T::Number) 
        i = ntuple(_->similar(λ) , 3)
        return  Dₜibb!(i , λ , T)
    end
        """
    Dₜibb(λ::AbstractVector,T::AbstractVector)

Calculates tuple of (`Ibb,dIbb/dT,d²Ibb/dT²`) calculated according to:

`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]`,
`d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]`,
`a₁=C₂/(λ*T)`,
`a₂ = 1/(eᵃ¹-1) ` , 
`a₃ = eᵃ¹/(eᵃ¹-1)` ,
`d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]`
    as far as 
        `Ibb = (λ⁻⁵)* C₁*a₂`
    and 
        `dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T)) = a₃*a₁*Ibb/T` 
    hense
        `d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
            = [a₃*a₁*Ibb/T^2]*[a₁*(2*a₃ - 1))-2] 
                = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2]` 
# Arguments:
`λ` - wavelength region, μm
`T` - temperature, Kelvins

# Returns:
`(Ibb , dIbb/dT , d²Ibb/dT²)`
    """
    function Dₜibb(λ::AbstractVector , T::AbstractVector{Q}) where Q
            # returns spectral intencity and its first and second derivatives with respect to the temperature
            (N , M) = (length(λ), length(T))
            i = ntuple(_-> Matrix{Q}(undef , N , M) , 3)
            return Dₜibb!(i , λ , T )
        end
        """
        band_power(T;λₗ=0.0 , λᵣ=Inf , tol=1e-8)

Total bb with temperature T integral intensity (`∫ibb(λ , T)dλ`) within (in-band radiance), [W/(m²⋅sr)]
the spectral range `λₗ...λᵣ` (by default the range is 0...inf)
tol - tolerance of intehration

# Arguments:

`T` - temperature,Kelvins
(optional)
`λₗ` - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
`tol` - intergation tolerance
    """
    band_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8) = power(T)*∫ibbₗ(T; λₗ=λₗ , λᵣ=λᵣ , tol=tol)


    """
    ∇ₜband_power(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)

Total bb with temperature T integral intensity derivative 
within (in-band radiance), [W/(m²⋅sr⋅K)]
the spectral range λₗ...λᵣ (by default the range is 0...inf)
tol - tolerance of integration

# Arguments:
T - temperature,Kelvins
(optional)
λₗ - left wavelength boundary, μm
λᵣ - right wavelength boundary, μm
tol - intergation tolerance
"""
∇ₜband_power(T; λₗ=0.0 , λᵣ=Inf , tol=1e-8) = ∇ₜpower(T) * ∫ibbₗ(T; λₗ=λₗ , λᵣ=λᵣ , tol=tol) + power(T) * ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )

"""
    ∇ₜband_power(T , band_power_value; λₗ=0.0 , λᵣ=Inf)

Evaluates band_power derivative with respect to temperature with `band_power_value`
provided externally (this may be usefull if one already has calculated the band_power itself)

# Arguments:
`T` - temperature,Kelvins
`band_power_value` - band_power calculated elswhere

(optional)
`λ`ₗ - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
`tol` - intergation tolerance

"""
function ∇ₜband_power(T , band_power_value; λₗ=0.0 , λᵣ=Inf) 
    P = power(T) 
    return (∇ₜpower(T) * band_power_value/P) + (P * ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ ))
end

 """
    ∇²ₜband_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)

Band power second derivative with respect to temeprature 

# Arguments:
`T` - temperature,Kelvins

(optional)
`λ`ₗ - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
`tol` - intergation tolerance
"""
function ∇²ₜband_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)
        return  ∇²ₜpower(T) * ∫ibbₗ(T; λₗ=λₗ , λᵣ=λᵣ , tol=tol) +
                2* ∇ₜpower(T) * ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ ) + 
                power(T) * ∇²ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )
    end

    """
    Dₜband_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)

returns all derivatives 
"""
function Dₜband_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)

        P = power(T)
        Pt = ∇ₜpower(T)
        Ptt = ∇²ₜpower(T)

        I = ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ , tol=tol) # this is the most time consuming
        It = ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )
        Itt = ∇²ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )

        bp  = P * I
        bpt = Pt * I + P * It
        bptt = Ptt * I + 2.0 * Pt * It + P * Itt

        return (bp , bpt , bptt)
    end    
    """
    Dₜband_power(T , skip_second_derivative::Val{true} ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)

returns the value and its derivative
"""
function Dₜband_power(T , skip_second_derivative::Val{true} ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)

        P = power(T)
        Pt = ∇ₜpower(T)

        I = ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ , tol=tol) # this is the most time consuming
        It = ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )

        bp  = P * I
        bpt = Pt * I + P * It

        return (bp , bpt , nothing)
    end    

  """
    ∇²ₜband_power(T , band_power_value; λₗ=0.0 , λᵣ=Inf , tol=1e-8)


Band power second derivative with respect to temeprature with `band_power` 
value evaluated externally 

# Arguments:
`T` - temperature,Kelvins

(optional)
`λ`ₗ - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
`tol` - intergation tolerance
"""
function ∇²ₜband_power(T , band_power_value; λₗ=0.0 , λᵣ=Inf)
    P = power(T)
    Pt = ∇ₜpower(T)

    # ∇ₜ (power * ∫ibbₗ) = ∇ₜpower * ∫ibbₗ + power * ∇ₜ∫ibbₗ
    # ∇ₜband_power(T , band_power_value) - ∇ₜpower * ∫ibbₗ = power * ∇ₜ∫ibbₗ
    # ∇ₜ∫ibbₗ = ( ∇ₜband_power(T , band_power_value) - ∇ₜpower * band_power_value/P)/P

    return  (∇²ₜpower(T) * band_power_value/P) +
                2.0 * Pt * (( ∇ₜband_power(T , band_power_value ; λₗ=λₗ , λᵣ=λᵣ) - Pt * band_power_value/P )/P) + 
                P * ∇²ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )
    end      
    
    """
    ∫ibbₗ(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)

Relative (with respect to the integral power in the whole spectrum)
integral intensity of bb in the spectral range `λₗ...λᵣ` 
(by default the range is 0...inf)

# Arguments:

`T - temperature,Kelvins

(optional)

`λₗ` - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
`tol` - intergation tolerance
"""
function ∫ibbₗ(T::Q; λₗ=0.0 , λᵣ=Inf , tol=1e-8) where Q
        # calculates the integral of spectral intencity over the wavelength
        @assert λₗ!=λᵣ "Bounding wavelengths must  not be equal"
        if λₗ>λᵣ
            (λₗ,λᵣ) = (λᵣ,λₗ)
        end
        if ~isfinite(λᵣ)# the right boundary is infinite
            if λₗ==0.0
                return one(Q)
            else #integration from fixed wavelength to infinity
                return ∫ibbₗ(T)- ∫ibbₗ(T , λᵣ=λₗ) 
            end
        else# righ wavelength boundary is finite
            if λₗ==0.0# integration from zero to fixed wavelength
                n=1
                ϵ = tol * 100
                summation = zero(Q)
                a = C₂/(λᵣ*T)
                while  (ϵ > tol) && (n < 1e4) 
                    etan = a*n
                    ϵ = (exp(-etan)/n)*(etan*(etan*(etan + 3.0) + 6.0) + 6.0)/(n^3) # there was a mistake 
                    summation += ϵ
                    n += 1;
                    #@show ϵ
                end
                return 15*summation/(pi^4)
            else# both wavelength resions are limited 
                return ∫ibbₗ(T , λᵣ=λᵣ) - ∫ibbₗ(T , λᵣ=λₗ)
            end
        end
    end


        """
        ∇ₜ∫ibbₗ(T; λₗ=0.0, λᵣ=Inf)

Relative (with respect to the integral power in the whole spectrum)
integral intensity derivative (analytic) of bb intensity fraction 
in the spectral range `λₗ...λᵣ` (by default the range is 0...inf)

# Arguments:
`T` - temperature,Kelvins
(optional)
`λₗ` - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
"""
    function ∇ₜ∫ibbₗ(T; λₗ=0.0, λᵣ=Inf)
                @assert λₗ != λᵣ "Bounding wavelengths must not be equal"
                if λₗ > λᵣ
                    (λₗ, λᵣ) = (λᵣ, λₗ)
                end
                
                xᵣ =  (λᵣ == 0.0 || !isfinite(λᵣ)) ? 0.0 : C₂ / (λᵣ * T)
                xₗ  =  (λₗ == 0.0  || !isfinite(λₗ)) ? 0.0 : C₂ / (λₗ * T)

                termᵣ = xᵣ == 0.0 ? 0.0 : (xᵣ^4) / expm1(xᵣ)
                termₗ = xₗ == 0.0 ? 0.0 : (xₗ^4) / expm1(xₗ)

                return (15.0 / (pi^4 * T)) * (termᵣ - termₗ)
        end

        """
    ∇²ₜ∫ibbₗ(T; λₗ=0.0, λᵣ=Inf)

Relative (with respect to the integral power in the whole spectrum)
integral intensity derivative (analytic) of bb intensity fraction 
in the spectral range `λₗ...λᵣ` (by default the range is 0...inf)

# Arguments:
`T` - temperature,Kelvins
(optional)
`λₗ` - left wavelength boundary, μm
`λᵣ` - right wavelength boundary, μm
    """
function ∇²ₜ∫ibbₗ(T; λₗ=0.0, λᵣ=Inf )
    @assert λₗ != λᵣ "Bounding wavelengths must not be equal"
    if λₗ > λᵣ
        (λₗ, λᵣ) = (λᵣ, λₗ)
    end

    xₗ = (λₗ == 0.0 || !isfinite(λₗ)) ? 0.0 : C₂ / (λₗ * T)
    xᵣ = (λᵣ == 0.0 || !isfinite(λᵣ)) ? 0.0 : C₂ / (λᵣ * T)

    termₗ = 0.0
    if xₗ > 0.0
        expm1_xₗ = expm1(xₗ)
        if isfinite(expm1_xₗ)
            termₗ = (xₗ^4 / expm1_xₗ) * (xₗ * exp(xₗ) / expm1_xₗ - 5.0)
        end
    end
    termᵣ = 0.0
    if xᵣ > 0.0
        expm1_xᵣ = expm1(xᵣ)
        if isfinite(expm1_xᵣ)
            termᵣ = (xᵣ^4 / expm1_xᵣ) * (xᵣ * exp(xᵣ) / expm1_xᵣ - 5.0)
        end
    end
    return (15.0 / (pi^4 * T^2)) * (termᵣ - termₗ)
end
        """
        attenuated_band_power(T, τ::AbstractVector, λ::AbstractVector;  tol=1e-6)

Evaluates the integrated band power of BB with temperature `T` reaching the detector, 
accounting for the spectral transmittance profile `τ(λ)`.

# Arguments:

`T` - temeperature , K 
`τ` - transmittance
`λ` - wavelength , μm
    """
    function attenuated_band_power(T, τ::AbstractVector, λ::AbstractVector)
        @assert length(λ) == length(τ) "Vectors λ and τ must have the same length"    
        return  planck_weighted(τ , λ , T)
    end
"""
    units(f::Function)

returns units string of output quantity  return 
"""
    function units(f::Function)  error(DomainError(f,"This function is not supported")) end

    units(::typeof(ibb)) = "W/(m²⋅sr⋅μm)" 
    units(::typeof(∇ₜibb)) = "W/(m²⋅sr⋅μm⋅K)"
    units(::typeof(∇²ₜibb)) = "W/(m²⋅sr⋅μm⋅K²)"
    units(::typeof(∇²ₗibb)) = "W/m²⋅sr⋅μm³"
    units(::typeof(power)) = "W/(m²⋅sr)"
    units(::typeof(band_power)) = "W/(m²⋅sr)"
    units(::typeof(λₘ)) = "μm"
    units(::typeof(tₘ)) = "K"



        """
        second_order_polynomial_fit(x1,x2,x3,g1,g2,g3)

    Hardcoded second order polynomial lsqr fitting
    """
    @inline function second_order_polynomial_fit(x1 , x2 , x3 , g1 , g2 , g3)

            d = (x1 - x2)*(x1 - x3)*(x2 - x3)
            inv_d = 1/d

            a1  = g3*x1^2*x2 - g2*x1^2*x3 - g3*x1*x2^2 + g2*x1*x3^2 + g1*x2^2*x3 - g1*x2*x3^2
            a2 = - g1*x2^2 + g2*x1^2 + g1*x3^2 - g3*x1^2 - g2*x3^2 + g3*x2^2
            a3 =  g1*x2 - g2*x1 - g1*x3 + g3*x1 + g2*x3 - g3*x2

            return (a1 * inv_d , a2 * inv_d , a3 * inv_d)
    end

    @inline fourth_order_polynomial_eval(a1, a2, a3, a4, a5, x) = @evalpoly(x, a1, a2, a3, a4, a5)

        """
    spectral_ratio(λ1::Number, λ2::Number, T::Number)

Calculate the theoretical intensity ratio `R =  ( Ibb1/Ibb2)`
between two wavelengths `λ1` and `λ2` at temperature `T`

# Arguments
`λ1`- First wavelength (usually the shorter one), in μm.
`λ2`- Second wavelength (usually the longer one), in μm.
`T`- Absolute temperature, in K.

"""
spectral_ratio(λ1::Number, λ2::Number, T::Number) = ibb(λ1 , T) /ibb(λ2 , T)

"""
    spectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;   tol = 1e-6) where TL <: Number

The same as spectral_ratio, but now the band can be wide (not a single wavelength).
This function may be useful for two-color pyrometers in case when their working regions 
width cannot be ignored 

# Arguments
`λ1`- tuple of left and right wavelength of the first band, in μm.
`λ2`- tuple of left and right wavelength of the second band, in μm.
`T`- Absolute temperature, in K.

"""
function spectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;   tol = 1e-6) where TL <: Number
   return band_power(T, λₗ = λ1[1] , λᵣ = λ1[2] , tol = tol) / band_power(T, λₗ = λ2[1] , λᵣ = λ2[2] , tol = tol)
end
"""
    ∇ₜspectral_ratio(λ1::Number, λ2::Number, T::Number)

Calculate the theoretical intensity ratio derivative `dR/T =d/dT (  Ibb1/Ibb2 )`
between two wavelengths `λ1` and `λ2` at temperature `T`

# Arguments
`λ1`- First wavelength (usually the shorter one), in μm.
`λ2`- Second wavelength (usually the longer one), in μm.
`T`- Absolute temperature, in K.
"""
function ∇ₜspectral_ratio(λ1::Number, λ2::Number, T::Number)
    (i1, di1 , _) = Dₜibb(λ1 , T)
    (i2, di2 , _) = Dₜibb(λ2 , T)
    return _spectral_ratio_first_derivative(i1 , di1 , i2 , di2)
end

"""
    ∇ₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;   tol = 1e-6) where TL <: Number

    First derivative of two wide spectral band ratio

# Arguments
`λ1`- tuple of left and right wavelength of the first band, in μm.
`λ2`- tuple of left and right wavelength of the second band, in μm.
`T`- Absolute temperature, in K
"""
function ∇ₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;   tol = 1e-6) where TL <: Number
    (i1, di1 , _) = Dₜband_power(T ; λₗ = λ1[1] , λᵣ = λ1[2] , tol = tol)
    (i2, di2 , _) = Dₜband_power(T ; λₗ = λ2[1] , λᵣ = λ2[2] , tol = tol)
    return _spectral_ratio_first_derivative(i1 , di1 , i2 , di2)
end
"""
    ∇²ₜspectral_ratio(λ1::Number, λ2::Number, T::Number)

Spectral ratio second derivative

# Arguments
`λ1`- First wavelength (usually the shorter one), in μm.
`λ2`- Second wavelength (usually the longer one), in μm.
`T`- Absolute temperature, in K.
"""
function ∇²ₜspectral_ratio(λ1::Number, λ2::Number, T::Number)
    (i1, di1 , d2i1) = Dₜibb(λ1 , T)
    (i2, di2 , d2i2) = Dₜibb(λ2 , T)
    return _spectral_ratio_second_derivative(
                i1, di1, d2i1,
                i2, di2, d2i2
    )
end
"""
    ∇²ₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;  tol = 1e-6) where TL <: Number

Spectral band ratio second derivative

# Arguments
`λ1`- tuple of left and right wavelength of the first band, in μm.
`λ2`- tuple of left and right wavelength of the second band, in μm.
`T`- Absolute temperature, in K.
"""
function ∇²ₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;  tol = 1e-6) where TL <: Number
    (i1, di1 , d2i1) = Dₜband_power(T, λₗ = λ1[1] , λᵣ = λ1[2] , tol = tol)
    (i2, di2 , d2i2) = Dₜband_power(T, λₗ = λ2[1] , λᵣ = λ2[2] , tol = tol)
    return  _spectral_ratio_second_derivative(
                i1, di1, d2i1,
                i2, di2, d2i2
    )
end
"""
    Dₜspectral_ratio(λ1::Number, λ2::Number, T::Number)

All spectral ratio derivatives in one tuple

# Arguments
`λ1`- First wavelength (usually the shorter one), in μm.
`λ2`- Second wavelength (usually the longer one), in μm.
`T`- Absolute temperature, in K.

"""
function Dₜspectral_ratio(λ1::Number, λ2::Number, T::Number)
    (i1, di1 , d2i1) = Dₜibb(λ1 , T)
    (i2, di2 , d2i2) = Dₜibb(λ2 , T)
    return (    
                 i1/i2 ,
                 _spectral_ratio_first_derivative(i1 , di1 , i2 , di2), 
                 _spectral_ratio_second_derivative(i1, di1, d2i1,
                                                            i2, di2, d2i2)
            )
end
function Dₜspectral_ratio(λ1::Number, λ2::Number, T::Number , skip_second_derivative::Val{true})
    (i1, di1 , _) = Dₜibb(λ1 , T , skip_second_derivative)
    (i2, di2 , _) = Dₜibb(λ2 , T , skip_second_derivative)
    return (    
                i1/i2 ,
                _spectral_ratio_first_derivative(i1 , di1 , i2 , di2), 
               nothing
            )
end
"""
    Dₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;  tol = 1e-6) where TL <: Number

All spectral band ratio derivatives are in one tuple


# Arguments
`λ1`- tuple of left and right wavelength of the first band, in μm.
`λ2`- tuple of left and right wavelength of the second band, in μm.
`T`- Absolute temperature, in K.
"""
function Dₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;  tol = 1e-6) where TL <: Number

    (i1, di1 , d2i1) = Dₜband_power(T, λₗ = λ1[1] , λᵣ = λ1[2] ; tol = tol)
    (i2, di2 , d2i2) = Dₜband_power(T, λₗ = λ2[1] , λᵣ = λ2[2] ; tol = tol)

    return (    
                i1/i2 ,
                _spectral_ratio_first_derivative(i1 , di1 , i2 , di2) , 
                _spectral_ratio_second_derivative(i1, di1, d2i1,
                                                            i2, di2, d2i2)
        )
end
    """
        Dₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number , skip_second_derivative::Val{true};  tol = 1e-6) where TL <: Number

    Ignores the second derivative evaluation 
    """
    function Dₜspectral_band_ratio(λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number , skip_second_derivative::Val{true};   tol = 1e-6) where TL <: Number

        (i1, di1 , _) = Dₜband_power(T, λₗ = λ1[1] , λᵣ = λ1[2] , skip_second_derivative; tol = tol)
        (i2, di2 , _) = Dₜband_power(T, λₗ = λ2[1] , λᵣ = λ2[2] , skip_second_derivative; tol = tol)

        return (    
                    i1/i2 ,
                    _spectral_ratio_first_derivative(i1 , di1 , i2 , di2),
                    nothing
            )
    end
    @inline """
    _spectral_ratio_second_derivative(
                    I1::Number, dI1_dT::Number, d2I1_dT2::Number,
                    I2::Number, dI2_dT::Number, d2I2_dT2::Number
        )

"""
function _spectral_ratio_second_derivative(
                    I1::Number, dI1_dT::Number, d2I1_dT2::Number,
                    I2::Number, dI2_dT::Number, d2I2_dT2::Number
        )
        I2_sq = I2 * I2
        inv_I2_cub = 1.0 / (I2_sq * I2) 
        
        numerator = d2I1_dT2 * I2_sq - 
                    I1 * I2 * d2I2_dT2 - 
                    2.0 * dI1_dT * dI2_dT * I2 + 
                    2.0 * I1 * (dI2_dT * dI2_dT)
                    
        return numerator * inv_I2_cub
    end
    @inline function _spectral_ratio_first_derivative(f1 , df1 , f2 , df2)
        return (df1 * f2 - df2 * f1 )/f2^2
    end
    ∇ₗ(::typeof(ibb)) = ∇ₗibb
    ∇ₗ(::typeof(∇ₗibb)) = ∇²ₗibb
    ∇²ₗ(::typeof(ibb)) = ∇²ₗibb
    # generating operators 
    for f in (:ibb, :band_power, :power, :spectral_ratio , :spectral_band_ratio , :planck_weighted , :planck_weighted_ratio)
        f_orig  = Symbol(f)         #  :ibb
        f_deriv = Symbol("∇ₜ", f)    #  :∇ₜibb
        f_sec   = Symbol("∇²ₜ", f)   #  :∇²ₜibb
        @eval begin
            ∇ₜ(::typeof($f_orig)) = $f_deriv
            ∇²ₜ(::typeof($f_orig)) = $f_sec
            ∇ₜ(::typeof($f_deriv)) = $f_sec
        end
    end
    ∫ₗ(::typeof(ibb)) = power
    ∫ₗ(::typeof(∇ₜibb)) = ∇ₜpower
    ∫ₗ(::typeof(∇²ₜibb)) = ∇²ₜpower
    struct BandIntegrator{F , L1 ,  L2}
        l1::L1
        l2::L2
        BandIntegrator{F}(l1::L1  , l2::L2) where {F , L1 , L2} = new{F , L1 , L2}(l1 , l2)
    end
    struct WeightedIntegrator{F , LT , AT}
        f::F 
        l::LT
        a::AT
    end
    (w::WeightedIntegrator)(t::Number) = weighted(w.a , w.l , t , w.f , identity )
    (b::BandIntegrator{typeof(ibb)})(T::Number) =   band_power(T , λₗ = b.l1 , λᵣ = b.l2)

    ∇ₜ(b::BandIntegrator{typeof(ibb)}) = BandIntegrator{typeof(∇ₜibb)}(b.l1 , b.l2)
    ∇²ₜ(b::BandIntegrator{typeof(ibb)}) = BandIntegrator{typeof(∇²ₜibb)}(b.l1 , b.l2)
    ∇ₜ(b::BandIntegrator{typeof(∇ₜibb)}) = BandIntegrator{typeof(∇²ₜibb)}(b.l1 , b.l2)

    (b::BandIntegrator{typeof(∇ₜibb)})(T::Number) =  ∇ₜband_power(T , λₗ = b.l1 , λᵣ = b.l2)
    (b::BandIntegrator{typeof(∇²ₜibb)})(T::Number) =   ∇²ₜband_power(T , λₗ = b.l1 , λᵣ = b.l2)
    struct RatioIntegrator{B1, B2}
        num::B1
        den::B2
    end
   /(int1::BandIntegrator, int2::BandIntegrator) = RatioIntegrator(int1, int2)
   /(int1::WeightedIntegrator, int2::WeightedIntegrator) = RatioIntegrator(int1, int2)

    (r::RatioIntegrator)(T::Number) = r.num(T) / r.den(T)
    """
    ∫ₗ(f::F , λ₁::Number, λ₂::Number) where F <: $(Dfunctions)

Returns a callable object which returns the f function integral as a function of temperature
```julia
f = ∫ₗ(ibb , 2.3 , 4.5)
f(1273.5) # returns the value of Planck function integrated over 2.3 - 4.5 spectral range at temperature 1273,5 K
```
"""
    ∫ₗ(::F, λ₁::Number, λ₂::Number) where F <: Dfunctions = BandIntegrator{F}(λ₁ , λ₂)
    ∫ₗ(g::Dfunctions , λ::AbstractVector , α::AbstractVector) =WeightedIntegrator( g , λ , α )
    """
        ∇ₜ(f) 

        Differentiation operator returns the derivatives function for the input function 
    [`ibb`](@ref) , [`band_power`](@ref) , [`power`](@ref) ,
    [`spectral_ratio`](@ref) , [`spectral_band_ratio`](@ref) , 
    [`planck_weighted`](@ref),[`planck_weighted_ratio`](@ref)
    and their derivatives e.g [`∇ₜband_power`](@ref) etc.

    # Examples
    ```julia
    ∇ₜ(ibb)        # -> ∇ₜibb
    F = (∇ₜ∘∇ₜ)(spectral_band_ratio)   #  ->  ∇²ₜspectral_band_ratio
    F(1253.5) # returns the second derivative with respect to temperature of the ratio of two spectral band at temperature 1253.5 K
    ```
    """
    function ∇ₜ end
      """
        ∇²ₜ(f) 

        Differentiation operator returns the derivatives of function
    [`ibb`](@ref) , [`band_power`](@ref) , [`power`](@ref) ,
    [`spectral_ratio`](@ref) , [`spectral_band_ratio`](@ref)

    # Examples
    ```julia
     ∇²ₜ(ibb)          # -> ∇²ₜibb
     ∇²ₜ(band_power)   # ->  ∇²ₜband_power
    ```
    """  
    function ∇²ₜ end        
    """
    ∫ₗ(g)
    ∫ₗ(g, λ₁, λ₂)
    ∫ₗ(g, λ::AbstractVector, α::AbstractVector)

The spectral integration operator for Planck radiation quantities.

Depending on the provided argument signature, this operator returns different 
highly optimized, type-stable callable objects (functors) designed to minimize 
runtime allocations.

# Arguments
* `g::Dfunctions`: The base Planck function layer (`ibb`, `∇ₜibb`, `∇²ₜibb`, or `Dₜibb`).
* `λ₁`, `λ₂::Number`: Spectral band boundaries in microns (μm).
* `λ::AbstractVector`, `α::AbstractVector`: Discrete experimental data arrays (wavelength grid and spectral coefficients).

---
### Full Spectrum Integration (Unary)
When applied directly to a single argument, it returns a functor representing the analytical total integrated power over the entire spectrum ([0, ∞)) governed by the Stefan-Boltzmann law.
```julia
f = ∫ₗ ibb
f(1273.5) # Returns σ * T^4
```

---

### Spectral Band Integration
When given a function and two numeric boundaries, it returns a functor that evaluates the integrated band power within the specific interval `[λ₁, λ₂]`.
```julia
f = ∫ₗ(ibb, 2.3, 4.5)
f(1273.5) # Returns integrated power between 2.3 and 4.5 μm
```

---
### Tabular Data Integration
When given a function and two matching discrete arrays, it returns a functor designed for experimental data pyrometry (e.g., integrating a wavelength-dependent emissivity vector `α` against the Planck distribution).
```julia
f = ∫ₗ(ibb, λ_grid, e_data)
f(1273.5) # Returns the discrete Planck-weighted integral 
```
"""
function ∫ₗ end        
    include("_chain_rules.jl")
    function symbolize end
    export symbolize

end


# Horner's scheme (explicit)
#=s1  += h * (c2_1 + h * (c3_1 + h * (c4_1 + h * c5_1)))
sn1 += h * (a1_1 + h * (a2_1 * inv2 + h * (a3_1 * inv3)))

s2  += h * (c2_2 + h * (c3_2 + h * (c4_2 + h * c5_2)))
sn2 += h * (a1_2 + h * (a2_2 * inv2 + h * (a3_2 * inv3)))

s3  += h * (c2_3 + h * (c3_3 + h * (c4_3 + h * c5_3)))
sn3 += h * (a1_3 + h * (a2_3 * inv2 + h * (a3_3 * inv3)))=#