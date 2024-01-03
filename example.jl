using Turing

@model function model(measurements_t, measurements_y)
    mass_quotient ~ Uniform(0.1, 10)
    observer_angle ~ Uniform(0., π/2)
    initial_phase ~ Uniform(-π, π)
    offset ~ Flat()

    predicted_magnitudes = star_magnitude(
        measurements_t, mass_quotient, observer_angle, ...
    )

    measurements_y .~ Normal.(predicted_magnitudes, σ)
end

m = model(measurements_t, measurements_y)
sample(m, 100000, NUTS())






