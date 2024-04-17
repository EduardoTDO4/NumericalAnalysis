function secant(f::Function, p_0::Float64, p_1::Float64, τ::Float64, N::Int64)

    #Step 1
    if f(p_0)==f(p_1)
        return("Erro, pois f(p_0) está muito próximo de f(p_1)")
    end

    #Step 2
    p = p_1; p_i = p_1; q_i = p_0; fp_i = f(p_i); fq_i = f(q_i)

    #Step 3
    for i = 1:N

        #Step 4
        if fq_i == fp_i
            return("Erro, pois f(p_i) está muito próximo de f(q_i)")
        end

        #Step 5
        p = p_i - (fp_i*(p_i - q_i))/(fp_i - fq_i)

        #Step 6
        if abs(p_i - q_i)<τ
            return (p)
        end

        #Step 7
        q_i = p_i; p_i = p; fp_i = f(p_i); fq_i = f(q_i)
    end
    
end