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

    [`∇²ₗibb`](@ref)  - spectral intensity second derivative with respect to wavelength

    [`rosseland_averaged_attenuation`](@ref) - Rosseland-averaging of spectral atternuation

    [`rosseland_averaged_attenuation`](@ref) - Planck-averaging of spectral atternuation

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

    Return tuple of three values:

    1 - a1 = C₂/(λ*T)
    2 - a2 = 1/expm1(a1)
    3 - a3 = exp(a1)*a2
"""
function  a₁₂₃(λ,T)
    a1=C₂/(λ*T)
    a2 = 1/expm1(a1)#1/eaxpm1(a)
    a3 = exp(a1)*a2#exp(a)/expm1(a)
    return isnan(a3) ? (a1,a2,0.0) : (a1,a2,a3)
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
function a₁₂₃!(amat::AbstractMatrix,λ::AbstractVector,T::Float64)
    # TRY VIEW WITH broadcating
        a1,a2,a3   = @views eachcol(amat)
        a1 .= (C₂/T)./λ
        a3 .=exp.(a1)
        a2 .= 1.0 ./(-1.0 .+a3)
        #a2 .= 1 ./expm1.(a1)
        a3 .= a3.*a2
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
    function ibb(λ,T) # general version close to the symbolic
        return (C₁/expm1(C₂/(λ*T)))*λ^-5
    end
    function ibb(λ::AbstractVector,T) # this version is useful in Optimization with AutoDiff?
        return map((l)->(C₁/expm1(C₂/(l*T)))*l^-5,λ)
    end
    function ibb(λ::AbstractVector,T::Base.RefValue{Float64}) # temeperature is a reference 
        return map((l)->(C₁/expm1(C₂/(l*T[])))*l^-5,λ)
    end
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
    function ibb(λ::AbstractVector,amat::AbstractMatrix) # internal version with provided coefficients matrix
        a2 = view(amat,:,2)
        return C₁*a2.*((1 ./λ).^5)     
    end
    """
    ibb(λ::AbstractVector,T::AbstractVector)

    Blackbody spectral intensity (spectral radiance),  [W/m²⋅sr⋅μm]
    Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1) , where a₁=C₂/(λ*T)
    Input:
        λ - wavelength in μm, [Nx0]
        T - temperature in Kelvins [Mx0]

"""
function ibb(λ::AbstractVector,T::AbstractVector)
        return length(T)==1 ? ibb(λ,T[1]) : @. ($C₁/expm1($C₂/(λ*$transpose(T))))*λ^-5
    end
    
    """
    ibb!(i::AbstractVector,λ::AbstractVector,T::Float64)
    
    In-place blackbody intensity,  [W/m²⋅sr⋅μm]
    Ibb = (λ⁻⁵)* C₁/(eᵃ¹-1) , where a₁=C₂/(λ*T)
    Input:
        i - bb intensity vector, [Nx0]
        λ - wavelength in μm, [Nx0]
        T - temperature in Kelvins
"""
function ibb!(i::AbstractVector,λ::AbstractVector,T::Float64)
        map!(l->(C₁/expm1(C₂/(l*T)))*l^-5,i,λ)
        return i
    end
    
"""
    ibb!(i::AbstractVector,λ::AbstractVector,amat::AbstractMatrix)::Nothing
    
    In-place blackbody intensity with intermediate coefficients provided externally, [W/m2-sr-mkm]
    Ibb =  C₁*(λ⁻⁵)*a₂ , where
    a₁=C₂/(λ*T)  - amat first column
    a₂ = 1/(eᵃ¹-1)  - amat second column 

    Input:
        i - BB intensity, [Nx0]
        λ - wavelength in μm,  [Nx0]
        amat - matrix of intermediate coefficients,  [Nx3]
"""
function ibb!(i::AbstractVector,λ::AbstractVector,amat::AbstractMatrix)::Nothing # this version is used in emissivity approximation 
        a2 = view(amat,:,2) # Ibb = (λ⁻⁵)* C₁*a₂
        i.=C₁*a2.*((1 ./λ).^5) 
        return  nothing   
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
    function ∇ₜibb(λ,T)
        return prod(a₁₂₃(λ,T))*C₁/(T*λ^5)
    end

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
function ∇ₜibb(λ::AbstractVector,T,amat::AbstractMatrix)
        return C₁*prod(amat;dims=2)./(T*λ.^5)
    end
    """
    ∇ₜibb!(g::AbstractMatrix,λ::AbstractVector,T::AbstractVector)

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
function ∇ₜibb!(g::AbstractMatrix,λ::AbstractVector,T::AbstractVector)
        # instance version of Planck function first derivative with respect to T
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                g[iii,jjj] = ∇ₜibb(l,t)
            end
        end   
    end
    function ∇ₜibb!(g::AbstractVector,λ::AbstractVector,T)# for fixed value of temperature
        # instance version of Planck function first derivative with respect to T
        #a = zeros(3)
        for (iii,l) in enumerate(λ) 
            g[iii] = prod(a₁₂₃(l,T))*C₁/(T*l^5)
        end 
    end
    # this version uses Vector for T because in this way the handle to the optimization varibale is implemented
    # T should be one-element array!
    """
    ∇ₜibb!(g::AbstractVector,λ::AbstractVector,T,amat::AbstractMatrix)

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
function ∇ₜibb!(g::AbstractVector,λ::AbstractVector,T,amat::AbstractMatrix)# for fixed value of temperature
        # instance version of Planck function first derivative with respect to T
        #a = zeros(3)
        #t = T[1]
        prod!(g,amat) # this puts all amat rows with column-wise product in g
        g./=(T*λ.^5) # a₃*a₁*a₂*C₁*(1/λ⁵)*(1/T)
        g.*=C₁
        return nothing
    end
    """
    ∇ₜibb!(g::AbstractVector,T, amat::AbstractMatrix,i::AbstractVector)::Nothing

    In-place bb intensity first derivative with respect to temperature
    with externally provided amat  - matrix with columns a₁,a₂,a₃

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
function ∇ₜibb!(g::AbstractVector,T, amat::AbstractMatrix,i::AbstractVector)::Nothing
        #t = T[1]
        a1 = view(amat,:,1)
        a3 = view(amat,:,3)
        g .= a1.*a3.*i/T  # dIbb/dT = a₃*a₁*Ibb/T
        return nothing
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
        a=a₁₂₃(λ,T)#ibb = C₁*a[2]*(l^-5)
        (a[1]*(2a[3]-1.0)-2.0)*a[1]*a[2]*a[3]*C₁/((T^3)*λ^5)
    end

    """
    ∇²ₜibb!(h::AbstractVector{Float64},λ::AbstractVector{Float64},T::Float64)

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
function ∇²ₜibb!(h::AbstractVector{Float64},λ::AbstractVector{Float64},T::Float64)# secpnd derivative for the fixed value of temperature
        # instance version of Planck function second derivative with respect to T
        for (iii,l) in enumerate(λ) 
            a=a₁₂₃(λ,T) # ibb = C₁*a[2]*(l^-5)         
            h[iii] = (a[1]*(2.0*a[3]-1.0)-2.0)*a[1]*a[2]*a[3]*C₁/((T^2)*l^5)
            # a₂*a₃*a₁*[a₁*(2*a₃ - 1))-2]*(C₁/(λ⁵*T²))
        end 
        return h
    end
    """
    ∇²ₜibb!(h::AbstractMatrix{Float64},λ::AbstractVector{Float64},T::AbstractVector{Float64})

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
function ∇²ₜibb!(h::AbstractMatrix{Float64},λ::AbstractVector{Float64},T::AbstractVector{Float64})
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                a=a₁₂₃(λ,T) # ibb = C₁*a[2]*(l^-5)         
                h[iii,jjj] = (a[1]*(2a[3]-1.0)-2.0)*a[1]*a[2]*a[3]*C₁/((t^2)*l^5) 
                # d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
            end
        end      
    end
    """
    ∇²ₜibb!(h::AbstractVector{Float64},λ::AbstractVector{Float64},T::Float64,amat::AbstractMatrix{Float64})::Nothing

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
function ∇²ₜibb!(h::AbstractVector{Float64},λ::AbstractVector{Float64},T::Float64,amat::AbstractMatrix{Float64})::Nothing
        # instance version of Planck function second derivative with respect to T
        # with supplied coefficints matrix
            a1 = view(amat,:,1)
            a3 = view(amat,:,3)
            prod!(h,amat) # h = a₃*a₁*a₂
            h./=(T^2)*λ.^5 # h = a₃*a₁*a₂*(1/λ⁵)*(1/T²)
            h.*=C₁ # h = C₁*a₃*a₁*a₂*(1/λ⁵)*(1/T²)
            h .*= a1.*(2.0*a3.-1) .-2.0  # C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
        return nothing        
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
function ∇²ₜibb!(h::AbstractVector{Float64},T::Float64,amat::AbstractMatrix{Float64},∇i::AbstractVector{Float64})::Nothing
        # instance version of Planck function second derivative with respect to T
        # with supplied coefficients matrix
            a1 = view(amat,:,1)
            a3 = view(amat,:,3)
            h.=∇i/T # h = (∇Ibb/T)
            h .*= a1.*(2.0*a3.-1.0) .-2.0  # h = (∇Ibb/T)*[a₁*(2*a₃ - 1))-2] 
        return nothing       
    end
    """
    rosseland_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number)

    Evaluates the Rosseland-averaged spectral attenuation coefficient (the summation of
    spectral scattering and absorption coefficients) α(λ) for temperature T:

    αᵣ = (∫(1/α(λ))⋅∇ₜibb(λ,T)dλ)⁻¹  
"""
rosseland_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number) = weighted_average(α,λ,T,∇ₜibb,inv)

"""
    planck_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number)

    Evaluates the Planck-averaged spectral attenuation coefficient (the summation of
    spectral scattering and absorption coefficients) α(λ) for temperature T:

    αᵣ = (∫(1/α(λ))⋅ibb(λ,T)dλ)⁻¹  
"""
planck_averaged_attenuation(α::AbstractVector, λ::AbstractVector,T::Number) = weighted_average(α,λ,T,ibb,inv) #identity

"""
    planck_averaged(x::AbstractVector, λ::AbstractVector,T::Number)

    Evaluates the Planck-averaged value of x(λ) for temperature T:

    xᵣ = ∫x(λ)ibb(λ,T)dλ 

    Can be used for example to evaluate the integral emissiovity from measured spectral emissivity
"""
planck_averaged(x::AbstractVector, λ::AbstractVector,T::Number) = weighted_average(x,λ,T,ibb,identity)
"""
    weighted_average(α::AbstractVector, λ::AbstractVector,T,g::Union{typeof(ibb),typeof(∇ₜibb),typeof(∇²ₜibb)},f::Function=identity)

    Generic function to evaluate the averaged value of some `fx` function of variable `x` dependent on wavelength `λ` 
    for temperature `T`. Uses linear approximation for the discrete variable and square polynomial 
    for the g function 

    xᵣ = f(∫f(x)g(λ,T)dλ), the default value of f is inv, e.g. if f = inv:
    xᵣ = 1/(∫g(λ,T)/x(λ)dλ)
"""
function weighted_average(α::AbstractVector, λ::AbstractVector,T,g::Union{typeof(ibb),typeof(∇ₜibb),typeof(∇²ₜibb)},f::Function=identity)
    @assert length(λ)==length(α)
    s = 0.0
    sn = 0.0
    norm = maximum(l->g(l,T),λ) # normalizing the value of Planck function
    g_T = l->g(l,T)/norm
    for i in eachindex(α)[begin:end-1]
        lstart = λ[i]
        lend = λ[i+1]
        lcentre = (lstart + lend )/2
        b1 = f(α[i]) 
        b2 = (f(α[i+1])-b1)/(lend-lstart)
        b1 = b1-b2*lstart
        #second order polynomial fitting of f within the interval
        (a1,a2,a3) = second_order_polynomial_fit(lstart,lcentre,lend,
                                                g_T(lstart),
                                                g_T(lcentre),
                                                g_T(lend))
        # g(l) = a1 + a2*l + a3*l²
        # f(l) = b1 + b2*l
        # ∫f(x)⋅g(x)dx = c2*l + c3*l² + c4*l³ + c5*l⁴
        c2=a1*b1
        c3=(a1*b2 + a2*b1)/2
        c4 = (a2*b2 + a3*b1)/3
        c5=a3*b2/4
        #evaluating the integrand values
        s+= fourth_order_polynomial_eval(0.0,c2,c3,c4,c5,lend) -
                            fourth_order_polynomial_eval(0.0,c2,c3,c4,c5,lstart)
        sn+=fourth_order_polynomial_eval(0.0,a1,a2/2,a3/3,0.0,lend) -
                         fourth_order_polynomial_eval(0.0,a1,a2/2,a3/3,0.0,lstart)
    end 
    return f(s/sn)
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
    ∇ₗibb(λ,T)


    BB intensity first derivative with respect to the wavelength

    Input:
        λ - wavelength, μm
        T - temperature, K
"""
    function ∇ₗibb(λ,T)
        # first derivative of Planck function with respect to wavelength
        #double a = C2/(lam*T);
        a=a₁₂₃(λ,T)
        return (a[1]/λ)*(a[3]-5/λ)*(C₁*a[2]*(λ^-5)) #(C₁*a₁₂₃(λ,T)[2])*λ^-5
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
        a=a₁₂₃(λ,T)
        return C₁*a[2]*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12)+30.0)/(λ^7)
    end
    """
    Dₗibb(λ,T)

    Returns a three-element tuple of (1.bb intensity,2.its first and 3.second derivative 
    with respect to the wavelentgh)

    Input:
        λ - wavelength, μm
        T - temperature, K
"""
function Dₗibb(λ,T)
        # methods returns PLanck function and its derivatives with respect to the wavelength
        # output is a tuple with (Planckfunction, Its first derivative with respect to the wavelength, Its second derivative with respect to the wavelength)
        a = a₁₂₃(λ,T)
        return (
            C₁*a[2]*((1/λ)^5),  # Planck function
            (C₁*a[2])*λ^-5,(a[1]/λ)*(a[3]-5/λ)*(C₁*a[2])*λ^-5, # first derivative
            (C₁/(λ^7) )*a[2]*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12)+30.0) # second derivative
        )
    end
    function Dₗibb(λ::AbstractVector,T::AbstractVector)
        # returns spectral intensity and its first and second derivatives with respect to the wavelength
        i = fill(0.0,length(λ), length(T))
        d1i = fill(0.0,length(λ), length(T))
        d2i = fill(0.0,length(λ), length(T))
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                a = a₁₂₃(l,t)
                i[iii,jjj] = C₁*a[2]*(l^-5)
                d1i[iii,jjj] = (a[1]/l)*(a[3]-5/l)*i[iii,jjj]
                d2i[iii,jjj] = (i[iii,jjj]/(l^2))*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12.0)+30.0)
                # (C₁/(λ^7))*a[2]*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12)+30)
            end
        end
        return (i,d1i,d2i)
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
        input_tuple, [Nx0 vector or nothing,Nx0 vector or nothing, Nx0 vector or nothing]
        λ - wavelength, μm, [Nx0]
        T - temperature, K   
""" 
    function Dₜibb!(input_tuple::Tuple{AbstractVector,AbstractVector,AbstractVector}, λ::AbstractVector,T)
        for (iii,l) in enumerate(λ) 
            a = a₁₂₃(l,T) #this function mutates global variable
            input_tuple[1][iii] = C₁*a[2]*(l^-5)   
            input_tuple[2][iii] = a[1]*a[3]*input_tuple[1][iii]/T #a[1]*a[2]*a[3]*C₁/(T*l^5)
            input_tuple[3][iii] = (a[1]*(2a[3]-1.0) -2.0)*input_tuple[2][iii]/T# = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
        end
        return input_tuple
    end
    function Dₜibb!(input_tuple:: Tuple{Nothing,AbstractVector,AbstractVector}, λ::AbstractVector,T)
        for (iii,l) in enumerate(λ) 
            a = a₁₂₃(l,T) #this function mutates global variable
            input_tuple[2][iii] = a[1]*a[3]*C₁*a[2]*(l^-5)/T #a[1]*a[2]*a[3]*C₁/(T*l^5)
            input_tuple[3][iii] = (a[1]*(2a[3]-1.0) -2.0)*input_tuple[2][iii]/T# = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
        end
        return input_tuple
    end
    function Dₜibb!(input_tuple:: Tuple{Nothing,Nothing,AbstractVector}, λ::AbstractVector,T::Float64)
        ∇²ₜibb!(input_tuple[3],λ,T)
        return input_tuple
    end
    function Dₜibb!(input_tuple:: Tuple{AbstractVector,Nothing,Nothing}, λ::AbstractVector,T::Float64)
        ibb!(input_tuple[2],λ,T)
        return input_tuple
    end
    function Dₜibb!(input_tuple::Tuple{Nothing,AbstractVector,Nothing},λ::AbstractVector,T::Float64)
        ∇ₜibb!(input_tuple[2],λ,T)
        return input_tuple
    end
    function Dₜibb!(input_tuple::Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}}, λ::AbstractVector,T::AbstractVector)
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                a = a₁₂₃(l,t) 
                input_tuple[1][iii,jjj] = C₁*a[2]*(l^-5)
                input_tuple[2][iii,jjj] = a[1]*a[3]*input_tuple[1][iii,jjj]/t #a[1]*a[2]*a[3]*C₁/(T*l^5)
                input_tuple[3][iii,jjj] = (a[1]*(2a[3]-1.0) -2.0)*input_tuple[2][iii,jjj]/t#(a[1]*(2a[3]-1)-2)*a[1]*a[2]*a[3]*C₁/((T^3)*l^5)
            end
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
        (Ibb,dIbb/dT,d²Ibb/dT²)
"""
function Dₜibb(λ::AbstractVector,T::AbstractVector)
        # returns spectral intencity and its first and second derivatives with respect to the temperature
        i = fill(0.0,length(λ), length(T))
        d1i = fill(0.0,length(λ), length(T))
        d2i = fill(0.0,length(λ), length(T))
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                a = a₁₂₃(l,t)
                i[iii,jjj] = C₁*a[2]*(l^-5)
                d1i[iii,jjj] = a[1]*a[3]*i[iii,jjj]/t #a[1]*a[2]*a[3]*C₁/(T*l^5)
                d2i[iii,jjj] = (a[1]*(2a[3]-1.0) -2.0)*d1i[iii,jjj]/t#(a[1]*(2a[3]-1)-2)*a[1]*a[2]*a[3]*C₁/((T^3)*l^5)
            end
        end
        return (i,d1i,d2i)
    end
    """
    band_power(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)

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
function band_power(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)
        return power(T)*∫ibbₗ(T;λₗ=λₗ,λᵣ=λᵣ,tol=tol)
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
function ∫ibbₗ(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)
        # calculates the integral of spectral intencity over the wavelength
        @assert λₗ!=λᵣ "Bounding wavelengths must be not equal"
        if λₗ>λᵣ
            (λₗ,λᵣ) = (λᵣ,λₗ)
        end
        if ~isfinite(λᵣ)# the right boundary is infinite
            if λₗ==0.0
                return 1
            else #integration from fixed wavelength to infinity
                return ∫ibbₗ(T)- ∫ibbₗ(T,λᵣ=λₗ) 
            end
        else# righ wavelength boundary is finite
            if λₗ==0.0# integration from zero to fixed wavelength
                n=1
                ϵ=tol*100
                summation=0.0
                a = C₂/(λᵣ*T)
                while  (ϵ>tol)&&(n<1e4) 
                    etan = a*n
                    summation+=(exp(-etan)/n)*(etan*(etan*(etan + 3.0) + 6.0) + 6.0)/(n^3) # there was a mistake 
                    n+=1;
                end
                return 15*summation/(pi^4)
            else# both wavelength resions are limited 
                return ∫ibbₗ(T,λᵣ=λᵣ) - ∫ibbₗ(T,λᵣ=λₗ)
            end
        end
    end
    """
    units(f::Function)

    returns units string of output quantity  return 
"""
function units(f::Function)  error(DomainError(f,"This function is unsupported")) end

units(::typeof(ibb)) = "W/m²⋅sr⋅μm" 
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
function second_order_polynomial_fit(x1,x2,x3,g1,g2,g3)
        d = (x1 - x2)*(x1 - x3)*(x2 - x3);
        a1  = g3*x1^2*x2 - g2*x1^2*x3 - g3*x1*x2^2 + g2*x1*x3^2 + g1*x2^2*x3 - g1*x2*x3^2
        a2 = - g1*x2^2 + g2*x1^2 + g1*x3^2 - g3*x1^2 - g2*x3^2 + g3*x2^2
        a3 =  g1*x2 - g2*x1 - g1*x3 + g3*x1 + g2*x3 - g3*x2
        return (a1/d,a2/d,a3/d)
    end
    function fourth_order_polynomial_eval(a1,a2,a3,a4,a5,x)
        return a1 + a2*x + a3*x^2 + a4*x^3+ a5*x^4
    end
end