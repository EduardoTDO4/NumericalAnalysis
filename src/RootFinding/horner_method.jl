function horner(p::Vector, p_0::Float64)

    n=length(p)-1

    if n == 0
        return(p[1],0.0)
    end

    if n == 1
        return (p[1]+p[2]*p_0,p_0)
    end

    a = p[n+1]
    b = p[n+1]

    for j = n:2
        a = p[j] + p_0*a
        b = a + p_0*b
    end
    a = p[1]+p_0*a
    return[a,b]
    
end

function newton_horner(p::Vector, p_0::Float64, τ::Float64, N::Int64)

    if length(p) == 1
        if p[1] == 0.0
            return("Toto número real é raíz de p")
        else
            return("p não possui raízes")
        end
    end

    if length(p) == 2
        if p[2] != 0.0
            return(-p[1]/p[2])
        elseif p[1] == 0.0
            return("Toto número real é raíz de p")
        else
            return("p não possui raízes")
        end
    end

    p_i = p_0
    q = p_0
 
    for i = 1:N
        
        w = horner(p,p_i)

        if w[2] == 0
            if w[1] == 0
                return(p_i)
            else
                return ("Método falhou em ", i)
            end
        end
 
        q = p_i - w[1] / w[2]

        if abs(q - p_i) < τ
            return (q)
        end
        p_i = q
 
    end
 
    return ("O método falhou")

end