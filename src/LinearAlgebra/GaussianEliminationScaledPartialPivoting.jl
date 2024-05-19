function pivoting(matriz_A::Matrix{Float64},vetor_b::Vector{Float64})

    n = length(vetor_b)

    dimension_failure(matriz_A,n) && return

    A = [copy(matriz_A) copy(vetor_b)]

    for i=1:n-1
        s = zeros(Float64,n)
        M_i = zeros(Float64,n)

        for k=i:n
            
            s[k] += maximum(abs,A[k,:])

            if s == zeros(Float64,n)
                return("Sistema não admite solução única")
            end

            M_i[k] += abs(A[k,i])/s[k]

        end

        M = maximum(abs,M_i)
        P = zeros(Int64,n-i+1)
        q=i
        
        for q = i:n
            if abs(A[q,i])/s[q] == M
                P[q-i+1] += q
            else
                P[q-i+1] += n
            end
        end
        
        Q = minimum(P)

        if Q!=i
            N = copy(A[Q,:])
            A[Q,:] = A[i,:]
            A[i,:] = N
        end
        
        for j = i+1:n
            m_ji = A[j,i]/A[i,i]
            A[j,:] -= m_ji*A[i,:]
        end     
    end

    if A[n,n] == 0
        return("Não há solução única")
    end

    x = zeros(Float64,n)
    x[n] += A[n,n+1]/A[n,n]

    for i = n-1:-1:1
        S = 0
        for j = i+1:n
            S+=A[i,j]*x[j]
        end
        x[i] += (A[i,n+1] - S)/A[i,i]
    end

    return(x)
end