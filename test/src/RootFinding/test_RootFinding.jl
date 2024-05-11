@testset "Test Root Findind" begin

    tol = 1e-10
    N = 100
    n = 30
    τ = tol

    f(x) = cos(x)
    Df(x) = -sin(x)
    f0 = Float64(pi / 2)
    g(x) = sin(x)
    Dg(x) = cos(x)
    g0 = Float64(pi)
    w(x) = (x^3) - 1
    Dw(x) = 3 * (x^2)
    w0 = 1.0
    y(x) = x^2
    Dy(x) = 2 * x
    y0 = 0.0
    z(x) = log(x)
    Dz(x) = 1 / x
    z0 = 1.0

    # bissection

    x_calc = bissection(f, 0.0, Float64(pi), τ, N)
    @test norm(f0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = bissection(g, Float64(pi / 2), Float64(3 * pi / 2), τ, N)
    @test norm(g0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = bissection(w, 0.0, 2.0, τ, N)
    @test norm(w0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = bissection(y, -1.0, 1.0, τ, N)
    @test norm(y0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = bissection(z, 0.1, 2.0, τ, N)
    @test norm(z0 - x_calc, Inf) ≈ 0 atol = tol

    # Fixed Point (só funciona com funções que são contrações, logo apenas com f e g)

    x_calc = fixed_point(f, 1.0, τ, N)
    @test norm(f0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = fixed_point(g, 2.5, τ, N)
    @test norm(g0 - x_calc, Inf) ≈ 0 atol = tol

    #x_calc = fixed_point(w,0.7,τ,N)
    #@test norm(w0-x_calc, Inf) ≈ 0 atol = tol

    #x_calc = fixed_point(y,0.2,τ,N)
    #@test norm(y0-x_calc, Inf) ≈ 0 atol = tol

    #x_calc = fixed_point(z,0.7,τ,N)
    #@test norm(z0-x_calc, Inf) ≈ 0 atol = tol

    # Newton

    x_calc = newton(f, Df, 1.0, τ, N)
    @test norm(f0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton(g, Dg, 2.5, τ, N)
    @test norm(g0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton(w, Dw, 0.7, τ, N)
    @test norm(w0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton(y, Dy, 0.2, τ, N)
    @test norm(y0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton(z, Dz, 0.7, τ, N)
    @test norm(z0 - x_calc, Inf) ≈ 0 atol = tol

    # Secant

    x_calc = secant(f, 0.0, 1.0, τ, N)
    @test norm(f0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = secant(g, 3.0, 3.5, τ, N)
    @test norm(g0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = secant(w, 0.0, 2.0, τ, N)
    @test norm(w0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = secant(y, -0.5, 0.2, τ, N)
    @test norm(y0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = secant(z, 0.1, 2.0, τ, N)
    @test norm(z0 - x_calc, Inf) ≈ 0 atol = tol

    # aitken (x^2 não satisfaz as hipóteses para aplicar aitken)

    x_calc = aitken(f,0.5,τ,n)
    @test norm(f0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = aitken(g,2.5,τ,n)
    @test norm(g0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = aitken(w,0.7,τ,n)
    @test norm(w0 - x_calc, Inf) ≈ 0 atol = tol

    #x_calc = aitken(y,0.5,τ,n)
    #@test norm(y0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = aitken(z,0.7,τ,n)
    @test norm(z0 - x_calc, Inf) ≈ 0 atol = tol

    # steffensen (como x^2 não funciona para aitken, ele não funciona para steffensen)

    x_calc = steffensen(f,1.0,τ,n)
    @test norm(f0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = steffensen(g,2.5,τ,n)
    @test norm(g0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = steffensen(w,0.7,τ,n)
    @test norm(w0 - x_calc, Inf) ≈ 0 atol = tol

    #x_calc = steffensen(y,0.7,τ,n)
    #@test norm(y0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = steffensen(z,0.7,τ,n)
    @test norm(z0 - x_calc, Inf) ≈ 0 atol = tol

    # horner (usei o método de horner para melhorar o método de newton)

    p = [1.0,1.0]
    p0 = -1.0
    q = [1.0,2.0,1.0]
    q0 = -1.0
    r = [0.0,0.0,1.0]
    r0 = 0.0
    s = [-1.0,0.0,0.0,1.0]
    s0 = 1.0
    t = [0.0,1.0,3.0,0.0,1.0]
    t0 = 0.0

    x_calc = newton_horner(p, 0.0, τ, N)
    @test norm(p0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton_horner(q, -0.5, τ, N)
    @test norm(q0 - x_calc, Inf) ≈ 0 atol = 1.0e-8

    x_calc = newton_horner(r, 0.7, τ, N)
    @test norm(r0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton_horner(s, 0.2, τ, N)
    @test norm(s0 - x_calc, Inf) ≈ 0 atol = tol

    x_calc = newton_horner(t, 0.7, τ, N)
    @test norm(t0 - x_calc, Inf) ≈ 0 atol = tol

end
