function aitken(f::Function,x_0::Float64,τ::Float64,N::Int64)

    h(x)=f(x)+x
    ϵ = 1.0e-16

    for i = 1:N
        x_1 = h(x_0)
        x_2 = h(x_1)
        denominador = (x_2 - x_1) - (x_1 - x_0)

        if abs(denominador)<ϵ
            return("Denominador muito próximo de zero")
        end

        aitkenx = x_2 - ( (x_2 - x_1)^2 )/denominador

        if abs(aitkenx - x_2)<τ
            return aitkenx
        end

        x_0 = aitkenx
    end
    return("Método falhou")
end