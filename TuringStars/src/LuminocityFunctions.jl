module LuminocityFunctions

using QuadGK
using Interpolations
using StructTypes

export
    T_4,
    planck_formula,
    black_body_K_rectangle,
    black_body_K,
    _black_body_K,
    claret_darkening

T_4(T) = T^4


function planck_formula(λ, T)
    h = 6.62607004e-34
    c = 299792458
    k = 1.38064852e-23
    return 2h*c^2 / (λ^5 * (exp(h*c / (λ*k*T)) - 1))
end

function black_body_K_rectangle(T)
    λ = 2.2e-6
    Δλ = 0.4e-6
    return planck_formula(λ, T) * Δλ
end



function make_interpolation_of_planck_integral(λ, Δλ, temperature_nodes)
    function luminocity_at_spectral_band(T)
        return quadgk(
            λ -> planck_formula(λ, T),
            λ - Δλ/2, λ + Δλ/2
        )[1]
    end
    luminocity_nodes = luminocity_at_spectral_band.(temperature_nodes)
    return cubic_spline_interpolation(temperature_nodes, luminocity_nodes)
end

const _black_body_K = make_interpolation_of_planck_integral(
    2.2e-6,
    0.4e-6,
    0 : 100 : 100_000
)

# JSON3.write(_black_body_K) = "Extrapolation([...very long lists...])

function black_body_K(T)
    return _black_body_K(T)
end

# JSON3.write(black_body_K) = "black_body_K"



StructTypes.StructType(::typeof(T_4)) = StructTypes.StringType()
StructTypes.StructType(::typeof(black_body_K_rectangle)) = StructTypes.StringType()
StructTypes.StructType(::typeof(black_body_K)) = StructTypes.StringType()


function claret_darkening(cosine, a1, a2, a3, a4)
    μ = cosine
    1. - a1 * (1 - √μ) - a2 * (1 - μ) - a3 * (1 - μ^1.5) - a4 * (1 - μ^2)
end

one1(x) = 1.

StructTypes.StructType(::typeof(claret_darkening)) = StructTypes.StringType()
StructTypes.StructType(::typeof(one1)) = StructTypes.StringType()

end