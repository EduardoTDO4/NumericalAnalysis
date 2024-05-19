@testset "Test inverting matrices" begin
    
    tol = 1.0e-10

    A = zeros(Float64,4,4)
    for i =1:4 
        A[i,i] += 1
    end
    B = A

    B_calc = inverse(A)
     
    for i = 1:4
        for j = 1:4
            @test norm(B[i,j]-B_calc[i,j],Inf) ≈ 0 atol = tol
        end
    end

    A = [ 2.0 3.0 ; 1.0 1.0]
    B = [-1.0 3.0 ; 1.0 -2.0]

    B_calc = inverse(A)
    for i = 1:2
        for j = 1:2
            @test norm(B[i,j]-B_calc[i,j],Inf) ≈ 0 atol = tol
        end
    end

end