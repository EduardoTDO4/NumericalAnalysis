function fixed_point(g::Function, p_0::Float64, τ::Float64, N::Int64)

    #Step 1
    h(x) = g(x) + x

    p = p_0
    p_i = p_0

    #Step 2
    for i = 1:N

        #Step 3
        p = h(p_i)

        #Step 4
        if abs(p - p_i) < τ
            return (p)
        end

        #Step 5
        p_i = p
    end

    #Step 6
    return ("Método falhou")

end
