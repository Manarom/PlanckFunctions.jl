# PlanckFunctions and its derivatives all tempeatures should be in Kelvins, all wavelengths in μm
module PlanckFunctions
    export  Dₜibb!,
            ∫ibbₗ,
            Tₖ, 
            ibb,
            ibb!,
            ∇ₜibb,
            ∇ₜibb!, 
            ∇²ₜibb,
            ∇²ₜibb!,
            ∇²ₗibb,
            a₁₂₃!,
            power,
            band_power,
            ∇ₜpower,
            ∇ₜband_power,
            planck_averaged,
            planck_averaged_attenuation,
            rosseland_averaged_attenuation,
            a₁₂₃,
            ∇ₜ,
            ∇²ₜ,
            ∇ₗ,
            ∇²ₗ

    const citation = "J.R.Howell,M.P.Menguc,J.R.Howell,M.P.Menguc,K.Daun,R.Siegel. Thermal radiation heat transfer. Seventh edition. 2021 " 
    const citation2 = "Risch, T.K., User's Manual: Routines for Radiative Heat Transfer and Thermometry. NASA/TM-2016-219103. 2016, Edwards, California: Armstrong flight Research Center"
    """
    `PlanckFunctions` module provides a set of functions for evaluating 
    the Planck thermal emission spectrum intensity (spectral radiance) and 
    its derivatives with respect to wavelength and temperature. It also provides 
    function to evaluate the integral over the wavelength radiance and Rosseland-averaged 
    and Planck-averaged of spectral coefficients. 

    Wavelength units are microns and all temperatures should be in Kelvins. 
    
    Main functions are:
     [`ibb`](@ref)  - spectral intensity (spectral radiance) in [W/m²⋅sr⋅μm]

     [`∇ₜibb`](@ref)  - spectral intensity first derivative with respect to temperature

    [`∇²ₜibb`](@ref)  - spectral intensity second derivative with respect to temperature

    [`∇ₗibb`](@ref)  - spectral intensity first derivative with respect to wavelength

    [`∇²ₗibb`](@ref)  - spectral intensity second derivative with respect to wavelength

    [`rosseland_averaged_attenuation`](@ref) - Rosseland-averaging of spectral attenuation

    [`planck_averaged_attenuation`](@ref) - Planck-averaging of spectral attenuation



    Main literature sources are:

    $(citation)
    
    $(citation2)
    """
    PlanckFunctions
    
    const ħ = 1.054_571_817E-34::Float64 # J*s
    #const C₁   = 1.191043E8::Float64#(1.191043E8,"W*μm*/m²*sr"," ","Risch","2016"),
    
    """
        C₁ constant for Planck function multiplier in [W⋅μm⁴/(m²⋅sr)]
        source $(citation) see Appendix A.
        Multiplied by 2 with respect to the source.
    """
    const C₁   = 2*0.595522001e8 # TRHT
    #const C₂ = 14387.752::Float64#(14387.752,"μm**K"," ","Risch","2016"),
    """
        C₂ constant for Planck spectral intensity exponent in [μm⋅K]
        source $(citation) see Appendix A.
    """
    const C₂ = 14387.7688::Float64# TRHT
    """
    C₃ constant of Wien's displacement law [μm⋅K]
    source $(citation) see Appendix A.
"""    
    const C₃ = 2897.77::Float64#(2897.77,"μm*K"," ","Risch","2016"),
    """
    C₄ constant in equation for maximum blackbody intensity [W/(m²⋅μm⋅sr*K⁵)]
    source $(citation) see Appendix A.
"""    
    const C₄ = 4.09567E-12::Float64 #TRHT
    #const C₄ = 4.09567E-12::Float64#(4.09567E-12,"W/m^2*μm*sr*K^5"," ", "Risch","2016"),
    """
    Stefan-Boltzmann constant [W/(m²*K⁴)]
    source $(citation) see Appendix A.
"""
    const σ  = 5.670367E-8::Float64
    #const σ  = 5.670400E-8::Float64 #(5.670400E-8,"W/(m²*K⁴)"," ", "Risch","2016"), previous constants values 
    const Tₖ = 273.15::Float64 #(273.15,"K"," ", "Risch","2016");
"""
    a₁₂₃(λ::Float64,T::Float64)

    Returns tuple of three values:

    1 - a1 = C₂/(λ*T)
    2 - a2 = 1/expm1(a1)
    3 - a3 = exp(a1)*a2
"""
function  a₁₂₃(λ::Number , T::F) where F <: Number
    a1 = C₂/(λ*T)
    a2 = 1.0/expm1(a1) # 1/eaxpm1(a)
    a3 = exp(a1)*a2  # exp(a)/expm1(a)
    return isnan(a3) ? (a1 , a2 , one(F)) : (a1 , a2 , a3)
end
"""
    a₁₂₃!(amat::AbstractMatrix,λ::AbstractVector,T::Float64)

    In-place filling of the intermediate matrix
    a₁=C₂/(λ*T)  - amat first column
    a₂ = 1/(eᵃ¹-1)  - amat second column 
    a₃ = eᵃ¹/(eᵃ¹-1) - amat third column

    Input:
        amat - matrix of intermediate coefficients size [Nx3]
        λ - wavelength in μm,  [Nx0]
        T - temperature in Kelvins
"""
function a₁₂₃!(amat::AbstractMatrix , λ::AbstractVector , T::Number)
    # TRY VIEW WITH broadcating
        a1 , a2 , a3   = @views eachcol(amat)
        @. a1 = (C₂/T)/λ
        @. a3 =exp(a1)
        @. a2 = 1.0 /(-1.0 + a3)
        #a2 .= 1 ./expm1.(a1)
        @. a3 = a3 * a2
    return amat
end
# BLACKBODY INTENSITY
"""
    ibb(λ,T)

    Blackbody spectral intensity (spectral radiance), [W/m²⋅sr⋅μm]
    Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1) , where a₁=C₂/(λ*T)

    Input:
        λ - wavelength in μm
        T - temperature in Kelvins
"""
    ibb(λ , T) = ( C₁/expm1( C₂/(λ*T) ) )*λ^-5
    ibb(λ::AbstractVector  , T::Number) = ibb.(λ , T)
    ibb(λ::AbstractVector,T::Base.RefValue{Q}) where Q <: Number = ibb(λ , T[])

"""
    ibb(λ::AbstractVector,amat::AbstractMatrix)

    Blackbody spectral intensity (spectral radiance)  with intermediate matrix provided externally, [W/m²⋅sr⋅μm]
    Ibb =  C₁*(λ⁻⁵)*a₂ , where
    a₁=C₂/(λ*T)  - amat first column
    a₂ = 1/(eᵃ¹-1)  - amat second column 
    Input:
        amat - matrix of intermediate coefficients,  [Nx3]
        λ - wavelength in μm,  [Nx0]
"""
    function ibb(λ::AbstractVector , amat::AbstractMatrix) # internal version with provided coefficients matrix
        a2 = @view amat[: , 2]
        return @. C₁ * a2 * ((1/λ)^5)     
    end
        """
        ibb(λ::AbstractVector,T::AbstractVector)

        Blackbody spectral intensity (spectral radiance),  [W/m²⋅sr⋅μm]
        Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1) , where a₁=C₂/(λ*T)
        Input:
            λ - wavelength in μm, [Nx0]
            T - temperature in Kelvins [Mx0]

    """
    ibb(λ::AbstractVector , T::AbstractVector) =  @. (C₁/expm1(C₂/(λ*$transpose(T))))*λ^-5

    
            """
        ibb!(i::AbstractVector , λ::AbstractVector , T::Number)
        
        In-place blackbody intensity,  [W/m²⋅sr⋅μm]
        Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1) , where a₁=C₂/(λ*T)
        Input:
            i - bb intensity vector, [Nx0]
            λ - wavelength in μm, [Nx0]
            T - temperature in Kelvins
    """
    ibb!(i::AbstractVector , λ::AbstractVector , T::Number) =  @. i = ibb(λ , T)

    
    """
        ibb!(i::AbstractVector , λ::AbstractVector , amat::AbstractMatrix)
        
    In-place blackbody intensity with intermediate coefficients provided externally, [W/m2-sr-mkm]
    Ibb =  C₁*(λ⁻⁵)*a₂ , where
    a₁=C₂/(λ*T)  - amat first column
    a₂ = 1/(eᵃ¹-1)  - amat second column 

    Input:
        i - BB intensity, [Nx0]
        λ - wavelength in μm,  [Nx0]
        amat - matrix of intermediate coefficients,  [Nx3]
    """
    function ibb!(i::AbstractVector , λ::AbstractVector , amat::AbstractMatrix) # this version is used in emissivity approximation 
            a2 = @view amat[: , 2] # Ibb = (λ⁻⁵)* C₁*a₂
            return @. i = C₁ * a2 *((1/λ)^5) 
        end   

    """
        ∇ₜibb(λ,T)

        BB intensity first derivative with respect to temperature
        dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))
        a₁=C₂/(λ*T)
        a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
        a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
        dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))
        Input:
            λ - wavelength in μm
            T - temperature in Kelvins
    """
    ∇ₜibb(λ , T) = prod(a₁₂₃(λ , T)) * C₁/(T * λ^5)
    ∇ₜibb(λ::AbstractVector , T::AbstractVector) = @. ∇ₜibb(λ , $transpose(T))    

    """
    ∇ₜibb(λ::AbstractVector,T,amat::AbstractMatrix)

    BB intensity first derivative with respect to temperature
    with externally provided matrix of intermediate coefficients
    dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   
    a₃ = eᵃ¹/(eᵃ¹-1) 
    dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))
    Input:
        λ - wavelength in μm, [Nx0]
        T - temperature in Kelvins
        amat - matrix of intermediate coefficients, [Nx3]
"""
function ∇ₜibb(λ::AbstractVector , T::Number , amat::AbstractMatrix)
        return C₁*prod(amat; dims=2)./(T*λ.^5)
    end
    """
    ∇ₜibb!(g::AbstractMatrix , λ::AbstractVector , T::AbstractVector)

    In-place BB intensity first derivative with respect to temperature
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))
    Input:
        g - vector to be filled, [Nx0]
        λ - wavelength in μm, [Nx0]
        T - temperature in Kelvins     
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

    dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))

    Input:
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

    Input:
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

    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]

    Input :
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

    Input :
    h  - to be filled, [Nx0]
    λ - wavelength in μm, [Nx0]
    T - tmperature in Kelvins  
"""
∇²ₜibb!(h::AbstractVector , λ::AbstractVector , T::Number)  = @. h =  ∇²ₜibb(λ,T)# secpnd derivative for the fixed value of temperature

    """
    ∇²ₜibb!(h::AbstractMatrix{Float64} , λ::AbstractVector{Float64}, T::AbstractVector{Float64})

    In-place bb intensity second order derivative with respect to temperature

    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]

    Input :
        h  - to be filled, [Nx0]
        λ - wavelength in μm, [Nx0]
        T - tmperature in Kelvins   
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

    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]

    Input :
    h  - to be filled, [Nx0]
    λ - wavelength in μm, [Nx0]
    T - temperature in Kelvins
    amat - matrix of intermediate coefficients,  [Nx3]
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

    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]
        as far as 
            Ibb = (λ⁻⁵)* C₁*a₂
        and 
            dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T)) = a₃*a₁*Ibb/T 
        hense
            d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
                = [a₃*a₁*Ibb/T^2]*[a₁*(2*a₃ - 1))-2] 
                    = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
    Input :
            h  - to be filled, [Nx0]
            λ - wavelength in μm, [Nx0]
            amat - matrix of intermediate coefficients,  [Nx3]
            ∇i - vector of bb intensity first derivatives, [Nx0]
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
    spectral scattering and absorption coefficients) α(λ) for temperature T:

    αᵣ = (∫(1/α(λ))⋅∇ₜibb(λ,T)dλ/∫∇ₜibb(λ,T)dλ)⁻¹  
"""
rosseland_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number) = inv(weighted_average(α,λ,T,∇ₜibb,inv))

"""
    planck_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number)

    Evaluates the Planck-averaged spectral attenuation coefficient (the summation of
    spectral scattering and absorption coefficients) α(λ) for temperature T:

    αᵣ = (∫(1/α(λ))⋅ibb(λ,T)dλ/∫ibb(λ,T)dλ)⁻¹  
"""
planck_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number) = inv(weighted_average(α,λ,T,ibb,inv)) #identity

"""
    planck_averaged(x::AbstractVector, λ::AbstractVector,T::Number)

    Evaluates the Planck-averaged value of x(λ) for temperature T:

    xᵣ = ∫x(λ)ibb(λ,T)dλ/∫ibb(λ,T)dλ 

    Can be used for example to evaluate the integral emissiovity from measured spectral emissivity
"""
planck_averaged(x::AbstractVector, λ::AbstractVector,T::Number) = weighted_average(x,λ,T,ibb,identity)
"""
    weighted_average(α::AbstractVector, 
                        λ::AbstractVector,
                        T , 
                        g::Union{typeof(ibb),typeof(∇ₜibb),typeof(∇²ₜibb)},
                        f::F = identity) where F <: Function

    Generic function to evaluate the averaged value of some `f(x)` function of variable `x` dependent
    on wavelength `λ` for temperature `T`. 
    Uses linear approximation for the discrete variable and square polynomial 
    for the g function 

    xᵣ = ∫f(x)g(λ,T)dλ/∫g(λ,T)dλ, the default value of f is `identity`,but, e.g. if f = inv:
    xᵣ = ∫g(λ,T)/x(λ)dλ/∫g(λ,T)dλ

"""
function weighted_average(α::AbstractVector, 
                        λ::AbstractVector,
                        T::Number , 
                        g::Union{typeof(ibb),typeof(∇ₜibb),typeof(∇²ₜibb)},
                        f::F = identity) where F <: Function
    N = length(α)
    @assert length(λ)==N

    s = 0.0
    sn = 0.0

    norm_val = maximum(l->g(l,T) , λ) # normalizing the value of Planck function
    nrm = 1/norm_val

    g_start = g(λ[begin], T) * nrm
    
    @inbounds for i in 1:N-1

        lstart = λ[i]
        lend = λ[i + 1]

        lcentre = (lstart + lend )/2

        # linear fitting of α
        b1 = f( α[i] ) 
        b2 = (f(α[i+1]) - b1)/(lend - lstart)
        b1 = b1 - b2 * lstart

        #second order polynomial fitting of f within the interval
        g_end = g(lend , T) * nrm #g_T(lend)
        g_centre = g(lcentre , T) * nrm
        (a1 , a2 , a3) = second_order_polynomial_fit(
                                                    lstart , lcentre , lend ,
                                                    g_start,
                                                    g_centre,
                                                    g_end
                                                )

        g_start = g_end                                        
        # g(l) = a1 + a2*l + a3*l²
        # f(l) = b1 + b2*l
        # ∫f(x)⋅g(x)dx = c2*l + c3*l² + c4*l³ + c5*l⁴
        c2=a1*b1
        c3=(a1*b2 + a2*b1)/2
        c4 = (a2*b2 + a3*b1)/3
        c5=a3*b2/4
        #evaluating the integrand values
        s += fourth_order_polynomial_eval(0.0 , c2 , c3 , c4 , c5 , lend) -
                            fourth_order_polynomial_eval(0.0 , c2 , c3 , c4 , c5 , lstart)
        sn += fourth_order_polynomial_eval(0.0 , a1 , a2/2 , a3/3 , 0.0 , lend) -
                         fourth_order_polynomial_eval(0.0 , a1 , a2/2 , a3/3 , 0.0 , lstart)
    end 
    return s/sn
end
    """
    λₘ(T)

    The wavelength (in μm) of bb intensity maximum vs temperature T 
    argmax(Planck(T))  - Wien's displacement law

    Input:
        T - temperature in Kelvins
"""
    function λₘ(T)
        # maximum wavelength of BB intencity in μm at temperature T (in Kelvins)
        C₃./T
    end
    """
    tₘ(λ)

    The temperature of BB having maximum at wavelength λ in Kelvins
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
function bright_temperature(i, λ; ϵ=1.0)
        return C₂ / (λ * log1p(ϵ * C₁ / (i * λ^5)))
    end
"""
    ∇ₗibb(λ,T)


    BB intensity first derivative with respect to the wavelength

    Input:
        λ - wavelength, μm
        T - temperature, K
"""
    function ∇ₗibb(λ,T)
        (a1 , a2 , a3) = a₁₂₃(λ,T)
        i = C₁ * a2 * (λ^-5)
        return (i/λ) * (a1 * a3 - 5.0)#
    end
    """
    ∇²ₗibb(λ,T)

    BB intensity second derivative with respect to the wavelength

    Input:
        λ - wavelength, μm
        T - temperature, K
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

    Input:
        λ - wavelength, μm
        T - temperature, K
"""
function Dₗibb(λ , T)
        # methods returns PLanck function and its derivatives with respect to the wavelength
        # output is a tuple with (Planckfunction, Its first derivative with respect to the wavelength, Its second derivative with respect to the wavelength)
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

    Returns integral (over the wavelength) intensity of BB (radiance)  at temperature T

    Units: W/(m²⋅sr)

    Input:
        T - temperature, K
"""
function power(T)
        # integral intencity of BB at temperature T
        return σ*(T^4)/π
    end
    """
    ∇ₜpower(T)

    Integral  (over the wavelength) intensity 
    derivative of BB (radiance)  at temperature T

    Units: W/(m²⋅sr⋅K)

    Input:
        T - temperature, K
"""
    ∇ₜpower(T) = 4*σ*(T^3)/π

    """
    ∇²ₜpower(T)

    Integral  (over the wavelength) intensity 
    second derivative of BB (radiance)  at temperature T

    Units: W/(m²⋅sr⋅K)

    Input:
        T - temperature, K
"""
∇²ₜpower(T) = 12.0 * σ * (T^2) / pi
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
        Dₜibb!(input_tuple, λ::AbstractVector,T)

        In-place filling the tuple of (bb intensity, its first ,and second ) derivatives with 
        respect to temperature
        d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
        d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
        a₁=C₂/(λ*T)
        a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
        a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
        d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]
            as far as 
                Ibb = (λ⁻⁵)* C₁*a₂
            and 
                dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T)) = a₃*a₁*Ibb/T 
            hense
                d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
                    = [a₃*a₁*Ibb/T^2]*[a₁*(2*a₃ - 1))-2] 
                        = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
        Input:
            input_tuple, {Nx0 vector or nothing , Nx0 vector or nothing, Nx0 vector or nothing}
            λ - wavelength, μm, [Nx0]
            T - temperature, K   
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
        """
        Dₜibb(λ::AbstractVector,T::AbstractVector)

        Calculates tuple of (Ibb,dIbb/dT,d²Ibb/dT²) calculated according to:
        d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
        d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
        a₁=C₂/(λ*T)
        a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
        a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
        d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]
            as far as 
                Ibb = (λ⁻⁵)* C₁*a₂
            and 
                dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T)) = a₃*a₁*Ibb/T 
            hense
                d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
                    = [a₃*a₁*Ibb/T^2]*[a₁*(2*a₃ - 1))-2] 
                        = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
        Input:
            λ - wavelength region, μm
            T - temperature, Kelvins
        Returns:
            (Ibb , dIbb/dT , d²Ibb/dT²)
    """
    function Dₜibb(λ::AbstractVector , T::AbstractVector{Q}) where Q
            # returns spectral intencity and its first and second derivatives with respect to the temperature
            (N , M) = (length(λ), length(T))
            i = ntuple(_-> Matrix{Q}(undef , N , M) , 3)
            return Dₜibb!(i , λ , T )
        end
        """
        band_power(T;λₗ=0.0 , λᵣ=Inf , tol=1e-8)

    Total bb with temperature T integral intensity within (in-band radiance), [W/(m²⋅sr)]
    the spectral range λₗ...λᵣ (by default the range is 0...inf)
    tol - tolerance of intehration

    Input:
        T - temperature,Kelvins
        (optional)
        λₗ - left wavelength boundary, μm
        λᵣ - right wavelength boundary, μm
        tol - intergation tolerance
    """
    band_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8) = power(T)*∫ibbₗ(T; λₗ=λₗ , λᵣ=λᵣ , tol=tol)


    """
    ∇ₜband_power(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)

        Total bb with temperature T integral intensity derivative 
    within (in-band radiance), [W/(m²⋅sr⋅K)]
    the spectral range λₗ...λᵣ (by default the range is 0...inf)
    tol - tolerance of integration

    Input:
        T - temperature,Kelvins
        (optional)
        λₗ - left wavelength boundary, μm
        λᵣ - right wavelength boundary, μm
        tol - intergation tolerance
"""
function ∇ₜband_power(T; λₗ=0.0 , λᵣ=Inf , tol=1e-8)
        return ∇ₜpower(T) * ∫ibbₗ(T; λₗ=λₗ , λᵣ=λᵣ , tol=tol) + power(T) * ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )
    end


 function ∇²ₜband_power(T ; λₗ=0.0 , λᵣ=Inf , tol=1e-8)
        return  ∇²ₜpower(T) * ∫ibbₗ(T; λₗ=λₗ , λᵣ=λᵣ , tol=tol) +
                2* ∇ₜpower(T) * ∇ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ ) + 
                power(T) * ∇²ₜ∫ibbₗ(T ; λₗ=λₗ , λᵣ=λᵣ )
    end   
    """
    ∫ibbₗ(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)

    Relative (with respect to the integral power in the whole spectrum)
    integral intensity of bb in the spectral range λₗ...λᵣ (by default the range is 0...inf)

    Input:
        T - temperature,Kelvins
        (optional)
        λₗ - left wavelength boundary, μm
        λᵣ - right wavelength boundary, μm
        tol - intergation tolerance
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
    integral intensity derivative of bb in the spectral range λₗ...λᵣ (by default the range is 0...inf)

    Input:
        T - temperature,Kelvins
        (optional)
        λₗ - left wavelength boundary, μm
        λᵣ - right wavelength boundary, μm
        tol - intergation tolerance
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
        integral intensity second derivative of bb in the spectral range λₗ...λᵣ (by default the range is 0...inf)

        Input:
            T - temperature,Kelvins
            (optional)
            λₗ - left wavelength boundary, μm
            λᵣ - right wavelength boundary, μm
            tol - intergation tolerance
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
    Input:
        T - temeperature , K 
        τ - transmittance
        λ - wavelength , μm
    """
    function attenuated_band_power(T, τ::AbstractVector, λ::AbstractVector;  tol=1e-6)
        @assert length(λ) == length(τ) "Vectors λ and τ must have the same length"    
        τ_avg = weighted_average(τ , λ , T , ibb , identity)
        return  τ_avg * band_power(T; λₗ=λ[begin], λᵣ=λ[end], tol=tol)
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


    ∇ₜ(::typeof(ibb)) = ∇ₜibb
    ∇ₜ(::typeof(∇ₜibb)) = ∇²ₜibb
    ∇²ₜ(::typeof(ibb)) = ∇²ₜibb
    ∇ₗ(::typeof(ibb)) = ∇ₗibb
    ∇ₗ(::typeof(∇ₗibb)) = ∇²ₗibb
    ∇²ₗ(::typeof(ibb)) = ∇²ₗibb

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

    @inline function fourth_order_polynomial_eval(a1, a2, a3, a4, a5, x)
        return @evalpoly(x, a1, a2, a3, a4, a5)
    end
end