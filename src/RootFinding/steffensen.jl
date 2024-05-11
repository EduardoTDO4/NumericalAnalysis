function steffensen(f::Function, x_0::Float64, τ::Float64, N::Int64)

    h(x) = f(x)+x

    for i = 1:N
        x_1=h(x_0)
        x_2=h(x_1)
        steffensen = x_0-(x_1-x_0)^2/(x_2-2*x_1+x_0)

        if abs(steffensen - x_0) < τ
            return(steffensen)
        end

        x_0 = steffensen

    end

    return("Método Falhou")
    
end