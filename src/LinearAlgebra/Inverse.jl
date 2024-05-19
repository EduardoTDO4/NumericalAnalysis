using LinearAlgebra

function inverse(matriz_A::Matrix{Float64})

    if size(matriz_A)[1]!=size(matriz_A)[2]
        return("A matriz deve ser quadrada")
    end

    if det(matriz_A) == 0
        return("A matriz não é inversível")
    end

    n = size(matriz_A)[1]
    B = zeros(Float64,n,n)
    C = zeros(Float64,n,n)
    for i =1:n 
        C[i,i] += 1
    end

    for i = 1:n
        B[:,i]+=pivoting(matriz_A,C[:,i])
    end

    return(B)
    
end