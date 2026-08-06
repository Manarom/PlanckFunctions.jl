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
    rrule(::typeof(ibb), λ::Number, T::Number)

Custom reverse-mode automatic differentiation rule (pullback) for `ibb(λ, T)`.
Safely intercepts the incoming adjoint `Δ` from upper execution layers. 
If `Δ` carries a `NoTangent` or `ZeroTangent` flag, it completely bypasses 
the calculation of analytical derivatives to save CPU cycles.
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
        
        ∂λ = ∇ₗibb(λ, T) * Δ
        ∂T = ∇ₜibb(λ, T) * Δ
        
        return ∂self, ∂λ, ∂T
    end
    
    return val, ibb_pullback
end