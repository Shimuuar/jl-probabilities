module LuminocityFunctions

using QuadGK
using Interpolations
using StructTypes

export
    T_4,
    planck_formula,
    black_body_K_rectangle,
    black_body_K,
    black_body_J,
    phoenixK,
    phoenixJ,
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


const _black_body_J = make_interpolation_of_planck_integral(
    1.25e-6,
    0.33e-6,
    0 : 100 : 100_000
)
function black_body_J(T)
    return _black_body_J(T)
end



const Ks = [6.739908535e15, 7.58105049e15, 8.412310905e15, 9.294638535e15, 1.022663698e16, 1.125340772e16, 1.23811021e16, 1.36827058e16, 1.52044639e16, 1.69739684e16, 1.886492455e16, 2.081903135e16, 2.2791595e16, 2.48251612e16, 2.69425769e16, 2.91309175e16, 3.116530765e16, 3.31222869e16, 3.48820174e16, 3.637262955e16, 3.756591015e16, 3.8675241e16, 3.96973334e16, 4.07117817e16, 4.171607145e16, 4.272045405e16, 4.37440985e16]

const Js = [1.572794644e16, 1.813480465e16, 2.062201313e16, 2.33238865e16, 2.621093235e16, 2.94456663e16, 3.301753605e16, 3.70857061e16, 4.15128043e16, 4.62876894e16, 5.1173646e16, 5.61922181e16, 6.12771178e16, 6.66556481e16, 7.26109868e16, 7.94075055e16, 8.63923993e16, 9.41957407e16, 1.023625334e17, 1.10382462e17, 1.182723065e17, 1.264122876e17, 1.346543525e17, 1.43024151e17, 1.51355723e17, 1.59907088e17, 1.68566452e17]

const Ts = 2300 : 100 : 4900


const _phoenixK = cubic_spline_interpolation(Ts, Ks, extrapolation_bc = Flat())
const _phoenixJ = cubic_spline_interpolation(Ts, Js, extrapolation_bc = Flat())

phoenixK(T) = _phoenixK(T)
phoenixJ(T) = _phoenixJ(T)




StructTypes.StructType(::typeof(T_4)) = StructTypes.StringType()
StructTypes.StructType(::typeof(black_body_K_rectangle)) = StructTypes.StringType()
StructTypes.StructType(::typeof(black_body_K)) = StructTypes.StringType()
StructTypes.StructType(::typeof(black_body_J)) = StructTypes.StringType()
StructTypes.StructType(::typeof(phoenixJ)) = StructTypes.StringType()
StructTypes.StructType(::typeof(phoenixK)) = StructTypes.StringType()


function claret_darkening(cosine, a1, a2, a3, a4)
    μ = cosine
    sqrt_μ = √μ
    I = one(μ)
    I - a1 * (I - sqrt_μ) - a2 * (I - μ) - a3 * (I - sqrt_μ^3) - a4 * (I - μ^2)
end

one1(x) = 1.

StructTypes.StructType(::typeof(claret_darkening)) = StructTypes.StringType()
StructTypes.StructType(::typeof(one1)) = StructTypes.StringType()

end