function ir_absolute_error(value::Irrational, approximated_value::Float64)::Float64
    # Implement the absolute error function. Hint: use abs funciton
    return abs(approximated_value - value)
end

function ir_relative_error(value::Irrational, approximated_value::Float64)::Float64
    # Implement the relative error function.
    @assert ~isapprox(value, 0, atol = 1e-8) "Value must be non-zero"
    return abs((value - approximated_value) / value)
end
