module Roche

using Meshes
using Roots

export Ω_potential, LagrangePoint_X, LagrangePoint_X_textbook, Ω_critical, roche_r, TransformToRocheLobe

function Ω_potential(r; mass_quotient, point_on_unit_sphere::Point)
    λ, μ, ν = coordinates(point_on_unit_sphere)
    return 1/r + 
        mass_quotient * (1 / √(1 + r^2 - 2r*λ) - λ*r) +
        (1. + mass_quotient) / 2 * r^2 * (1. - ν^2)
end

function LagrangePoint_X(mass_quotient)
    return 1 / (1. + mass_quotient)
end

function LagrangePoint_X_textbook(mass_quotient)
    function Ω_x_derivative(x)
        -1 / x^2 - 
        mass_quotient * (1. - 1 / (x-1)^2) + 
        (1 + mass_quotient) * x
    end
    return find_zero(Ω_x_derivative, (0., 1.))
end

function Ω_critical(lagrange1_x, mass_quotient)
    return Ω_potential(lagrange1_x; mass_quotient, point_on_unit_sphere=Point(1., 0., 0.))
end

function roche_r(Ω0, lagrange1_x, mass_quotient, point_on_unit_sphere::Point)
    function Ω_partially_applied(r)
        return Ω_potential(r; mass_quotient, point_on_unit_sphere) - Ω0
    end
    find_zero(Ω_partially_applied, (0., lagrange1_x))
end


# Similar to https://github.com/JuliaGeometry/Meshes.jl/blob/a0487c6824d6ee9d7389edc25ae937f1e4cf26fd/src/transforms/translate.jl

struct TransformToRocheLobe <: StatelessGeometricTransform
    mass_quotient
    lagrange1_x
    Ω0
end

function TransformToRocheLobe(mass_quotient)
    lagrange1_x = LagrangePoint_X_textbook(mass_quotient)
    Ω0 = Ω_critical(lagrange1_x, mass_quotient)
    return TransformToRocheLobe(mass_quotient, lagrange1_x, Ω0)
end

Meshes.preprocess(transform::TransformToRocheLobe, object) = transform

function Meshes.applypoint(::TransformToRocheLobe, points, prep::TransformToRocheLobe)

    function transform_point(point::Point)
        r = roche_r(prep.Ω0, prep.lagrange1_x, prep.mass_quotient, point)
        Point(coordinates(point) .* r)
    end

    return map(transform_point, points), prep
end


end