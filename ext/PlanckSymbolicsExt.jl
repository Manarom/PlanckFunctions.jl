module PlanckSymbolicsExt

    import PlanckFunctions as PF
    using Symbolics 

    PF.symbolize(f::PF.TwoArgsFunctions , λ::Num, T::Num) = f(λ, T)
    PF.symbolize(f::PF.RatioFunctions , λ1::Num, λ2::Num ,  T::Num)  = f(λ1 , λ2 , T)
    PF.symbolize(f ,  T::Num) = begin 
            try
                return f(T)
            catch e2
                error("Symbolization is not possible for function $f. Error: $e2")
            end
        end

end