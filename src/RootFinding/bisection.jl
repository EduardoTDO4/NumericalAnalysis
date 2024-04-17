function bissection(f::Function,a::Float64,b::Float64,τ::Float64, N::Int64)

    #Step 1
    a_i=a; b_i=b

    #Step 2
    for i=1:N

        #Step 3
        p = a_i + (b_i - a_i)/2; f_p = f(p)

        #Step 4
        if f_p == 0 | (b_i - a_i)/2 < τ
            return(p);
        end

        #Step 5
        if sign(f(a_i))*sign(f_p)>0
            a_i = p
        else
            b_i = p
        end

    end

    #Step 6 
    return("Método falhou")
    
end