module QuadGKPlanckExt

using PlanckFunctions
using QuadGK
import PlanckFunctions: ∫ₗ, ContinuousTemperatureWeightedIntegrator, Dₜibb

# 1. High-Performance Combined Vector-Quadrature Evaluation
function Dₜplanck_weighted(q::AbstractSpectralQuantity, λ₁::Number, λ₂::Number, T::Number; rtol=1e-8, atol=1e-12)
    # The integrand function evaluates all 3 Leibniz layers simultaneously.
    # QuadGK integrates this tuple efficiently by adapting its grid over the combined error profile.
    function integrand(λ)
        ibb_val, d_ibb, d2_ibb = Dₜibb(λ, T)
        q_val, d_q, d2_q = PlanckFunctions.Dₜ(q, λ, T)
        
        # Leibniz integral rule expansions
        f_val  = q_val * ibb_val
        f_grad = q_val * d_ibb  +  d_q * ibb_val
        f_hess = q_val * d2_ibb + 2.0 * d_q * d_ibb + d2_q * ibb_val
        
        return (f_val, f_grad, f_hess)
    end

    integrals, _ = quadgk(integrand, λ₁, λ₂; rtol=rtol, atol=atol)
    return integrals
end

# 2. Overriding the Functor to activate the real QuadGK loop
function (w::ContinuousTemperatureWeightedIntegrator)(T::Number)
    return continuous_Dₜweighted_value(w.quantity, w.l1, w.l2, T)
end

# 3. Extending the DSL operator dispatch
# Direct interface implementation
PlanckFunctions.∫ₗ(q::AbstractSpectralQuantity, λ₁::Number, λ₂::Number) = 
    ContinuousTemperatureWeightedIntegrator(q, λ₁, λ₂)

# Implicit wrapper for standard user functions (automatically applies ForwardDiff layer)
PlanckFunctions.∫ₗ(ϵ_func::Function, λ₁::Number, λ₂::Number) = 
    ContinuousTemperatureWeightedIntegrator(AutodiffSpectralQuantity(ϵ_func), λ₁, λ₂)

end # module PlanckQuadGKExt