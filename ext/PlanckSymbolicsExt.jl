module PlanckSymbolicsExt

    import PlanckFunctions as PF
    using Symbolics 

    function PF.symbolize(f::Function)
        Symbolics.@variables λ T
    
        if hasmethod(f, Tuple{Any, Any})
            return f(λ, T)
        end
        
        if hasmethod(f, Tuple{Any})
            return f(T)
        end
        
        error("Symbolization is not possible for this function")

    end

end