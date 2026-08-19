module PlanckQuadGKExt

import PlanckFunctions as PF
using QuadGK , StaticArrays

function PF.Dₜplanck_weighted(q::PF.AbstractSpectralQuantity, λ₁::Number, λ₂::Number, T::Number; rtol=1e-8, atol=1e-12)

    function integrand(λ)
        ibb_val, d_ibb, d2_ibb = PF.Dₜibb(λ, T)
        q_val, d_q, d2_q = PF.eval_Dₜ(q, λ, T)
        
        # Leibniz integral rule expansions
        f_val  = q_val * ibb_val
        f_grad = q_val * d_ibb  +  d_q * ibb_val
        f_hess = q_val * d2_ibb + 2.0 * d_q * d_ibb + d2_q * ibb_val
        
        return SVector(f_val, f_grad, f_hess)
    end

    integrals, _ = quadgk(integrand, λ₁, λ₂; rtol=rtol, atol=atol)
    return Tuple(integrals)
end
function PF.planck_weighted(q::PF.AbstractSpectralQuantity,
         λ₁::Number, 
         λ₂::Number, 
         T::Number; rtol=1e-8, atol=1e-12)
    function integrand(λ)
        ibb_val = PF.ibb(λ, T)
        q_val = q(λ, T)
        f_val  = q_val * ibb_val
        return f_val
    end
    val, _ = quadgk(integrand, λ₁, λ₂; rtol=rtol, atol=atol)
    return val
end

function PF.∇ₜplanck_weighted(q::PF.AbstractSpectralQuantity, λ₁::Number, λ₂::Number, T::Number; rtol=1e-8, atol=1e-12)
    function integrand(λ)
        (ibb_val , d_ibb, _) = PF.Dₜibb(λ, T)
        (q_val , d_q , _) = PF.eval_Dₜ(q, λ, T)
        f_grad = q_val * d_ibb  +  d_q * ibb_val 
        return f_grad
    end
    val, _ = quadgk(integrand, λ₁, λ₂; rtol=rtol, atol=atol)
    return val
end
function PF.∇²ₜplanck_weighted(q::PF.AbstractSpectralQuantity, λ₁::Number, λ₂::Number, T::Number; rtol=1e-8, atol=1e-12)
    function integrand(λ)
        ibb_val, d_ibb, d2_ibb = PF.Dₜibb(λ, T)
        q_val, d_q, d2_q = PF.eval_Dₜ(q, λ, T)
        f_hess = q_val * d2_ibb + 2.0 * d_q * d_ibb + d2_q * ibb_val
        return f_hess
    end
    val, _ = quadgk(integrand, λ₁, λ₂; rtol=rtol, atol=atol)
    return val
end

end # module PlanckQuadGKExt