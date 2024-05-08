function ir_absolute_error(value::Irrational, approximated_value::Float64)::Float64
    abs(value - approximated_value)
end
