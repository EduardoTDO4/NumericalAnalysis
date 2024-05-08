function newton(f::Function, Df::Function, p_0::Float64, τ::Float64, N::Int64)

    #Step 1
    p_i = p_0
    p = p_0

    #Step 2
    for i = 1:N
        #Step 3
        if Df(p_i) == 0
            return ("Método falhou")
        end

        #Step 4
        p = p_i - f(p_i) / Df(p_i)

        #Step 5
        if abs(p - p_i) < τ
            return (p)
        end

        #Step 6
        p_i = p

    end

    #Step 7
    return ("O método falhou")

end
