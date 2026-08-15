module PlanckSymbolicsExt

    import PlanckFunctions as PF
    using Symbolics 

    function PF.symbolize(f::Function)
        # Создаем символьные переменные
        Symbolics.@variables λ T
    
        try
            return f(λ, T)
        catch e
            try
                return f(T)
            catch e2
                error("Symbolization is not possible for function $f. Error: $e2")
            end
        end
    end

end