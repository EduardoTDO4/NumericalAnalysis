@testset "Test inverting matrices" begin
    
    tol = 1.0e-10

    A = eye(4)
    B = eye(4)

    B_calc = inverse(A)
    @test 
    for i = 1:4
        for j = 1:4
            norm(B[i,j]-B_calc[i,j],inf) == 0 atol=tol
        end
    end

end