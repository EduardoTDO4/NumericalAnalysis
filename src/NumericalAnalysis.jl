module NumericalAnalysis

include("Includes.jl")

export round_sum,
    round_mul,
    trunc_sum,
    trunc_mul,
    absolute_error,
    relative_error,
    approximated_value_interval_given_relative_error,
    bissection,
    fixed_point,
    newton,
    secant,
    solve_system
    aitken
    steffensen
    horner
    newton_horner
    dimension_failure
    pivoting
    inverse
end
