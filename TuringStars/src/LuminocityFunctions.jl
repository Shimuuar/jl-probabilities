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



const Ks = [1.8703479891078204e16, 2.1263852966397804e16, 2.3823204641797948e16, 2.6533323297406332e16, 2.93597771009047e16, 3.242922313466197e16, 3.5745423654219784e16, 3.9493620073171464e16, 4.379088435379438e16, 4.877045094675458e16, 5.416349622314456e16, 5.980896723520685e16, 6.5491249472586104e16, 7.122249385435577e16, 7.70488245890636e16, 8.296050494690677e16, 8.848011316502115e16, 9.382427626724853e16, 9.873491505421731e16, 1.0297676274742496e17, 1.0642306717787766e17, 1.0965276642183946e17, 1.1273511198227566e17, 1.1579283551593573e17, 1.188218235158687e17, 1.2185226356103037e17, 1.249341260982901e17, 1.2786868272007986e17]

const Js = [4.22043792077667e16, 4.943081244273363e16, 5.6970008900466424e16, 6.515104885509999e16, 7.383626419624864e16, 8.347352683645163e16, 9.402297836551931e16, 1.0589221933904166e17, 1.1864196531573013e17, 1.3219732908270312e17, 1.4588751053141686e17, 1.597073720345217e17, 1.7335401619788986e17, 1.8727378627642288e17, 2.02082331754589e17, 2.183574817037801e17, 2.3496231409447325e17, 2.5353092083034022e17, 2.7326537448986157e17, 2.935292732184174e17, 3.1500129641145574e17, 3.374266401802471e17, 3.605229332900457e17, 3.8413868191222195e17, 4.077673628117195e17, 4.320874411492431e17, 4.568033760719095e17, 4.813481566833534e17]

const Ts = 2300 : 100 : 5000


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
    1. - a1 * (1 - sqrt_μ) - a2 * (1 - μ) - a3 * (1 - sqrt_μ^3) - a4 * (1 - μ^2)
end

one1(x) = 1.

StructTypes.StructType(::typeof(claret_darkening)) = StructTypes.StringType()
StructTypes.StructType(::typeof(one1)) = StructTypes.StringType()

end