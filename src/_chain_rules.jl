# chainrules for AD 
## FORWARD
"""
    frule((Δself, Δλ, ΔT), ::typeof(ibb), λ::Number, T::Number)

Custom forward-mode automatic differentiation rule for the Planck function `ibb(λ, T)`
"""
function ChainRulesCore.frule(
                                (_, Δλ, ΔT), 
                                ::typeof(ibb), 
                                λ::Tλ, T::TT
                            ) where {Tλ <:Number , TT<:Number}
    val = ibb(λ, T)
    ∂val = ZeroTangent()
    if !(Δλ isa NoTangent || Δλ isa ZeroTangent)
        Δλ_val = unthunk(Δλ) 
        ∂val += ∇ₗibb(λ, T) * Δλ_val
    end
    
    if !(ΔT isa NoTangent || ΔT isa ZeroTangent)
        ΔT_val = unthunk(ΔT) 
        ∂val += ∇ₜibb(λ, T) * ΔT_val
    end
    return val, ∂val
end
"""
    frule((Δself, ΔT), ::typeof(band_power), T; λₗ=0.0, λᵣ=Inf, tol=1e-8)

Forward-mode AD rule for `band_power`. 
"""
function ChainRulesCore.frule(
                        (_, ΔT), 
                        ::typeof(band_power), 
                        T::Number; 
                        λₗ=0.0, λᵣ=Inf, tol=1e-8
                    )
    val = band_power(T; λₗ=λₗ, λᵣ=λᵣ, tol=tol)
    ∂val = ZeroTangent()
    
    if !(ΔT isa NoTangent || ΔT isa ZeroTangent)
        ΔT_val = unthunk(ΔT)
        ∂val += ∇ₜband_power(T, val; λₗ=λₗ, λᵣ=λᵣ) * ΔT_val
    end
    
    return val, ∂val
end

"""
    frule((Δself, ΔT), ::typeof(band_power), T; λₗ=0.0, λᵣ=Inf, tol=1e-8)

Forward-mode AD rule for `band_power`. 
"""
function ChainRulesCore.frule( (_ , _ , _ , ΔT), 
                        ::typeof(spectral_ratio),
                        λ1::Number ,  λ2::Number ,
                        T::Number)

    (val , dR, _) = Dₜspectral_ratio(λ1 , λ2 , T , Val(true))

    ∂val = ZeroTangent()
    
    if !(ΔT isa NoTangent || ΔT isa ZeroTangent)
        ΔT_val = unthunk(ΔT)
        ∂val += dR * ΔT_val
    end
    
    return val, ∂val
end

# spectral_band_ratio  # (λ1::NTuple{2, TL}, λ2::NTuple{2,TL}, T::Number;  tol = 1e-6)
"""
    rrule(::typeof(ibb), λ::Number, T::Number)

Custom reverse-mode automatic differentiation rule (pullback) for `ibb(λ, T)`.
"""
function ChainRulesCore.rrule(
                        ::typeof(ibb), 
                        λ::Tλ, T::TT
                    ) where {Tλ<:Number, TT<:Number}
    val = ibb(λ, T)
    function ibb_pullback(Δ)
        ∂self = NoTangent() 
        
        if Δ isa NoTangent
            return (∂self, NoTangent(), NoTangent())
        elseif Δ isa ZeroTangent
            return (∂self, ZeroTangent(), ZeroTangent())
        end
        
        ∂λ = @thunk(∇ₗibb(λ, T) * Δ)
        ∂T = @thunk(∇ₜibb(λ, T) * Δ)
        
        return ∂self, ∂λ, ∂T
    end
    
    return val, ibb_pullback
end

"""
    ChainRulesCore.rrule(
                        ::typeof(band_power), 
                        T::TT; λₗ::Number , λᵣ::Number
                    ) where { TT<:Number}


Custom reverse-mode automatic differentiation rule (pullback) for `band_power(T; λₗ = λₗ , λᵣ = λᵣ)`.
"""
function ChainRulesCore.rrule(
                        ::typeof(band_power), 
                        T::TT; λₗ::Number , λᵣ::Number
                    ) where { TT<:Number}
    val = band_power(T ; λₗ = λₗ , λᵣ = λᵣ)
    function _pullback(Δ)
        ∂self = NoTangent() 
        
        if Δ isa NoTangent
            return (∂self, NoTangent())
        elseif Δ isa ZeroTangent
            return (∂self, ZeroTangent())
        end
        ∂T = @thunk(∇ₜband_power(T ; λₗ = λₗ , λᵣ = λᵣ) * Δ)
        
        return ∂self,  ∂T
    end
    
    return val, _pullback
end