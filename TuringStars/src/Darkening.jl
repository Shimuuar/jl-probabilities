module Darkening

using Interpolations
using DelimitedFiles
using StructTypes

export
    claret_coefs_interpolant,
    J_coefs_interpolant,
    K_coefs_interpolant

function claret_coefs_interpolant(channel::Symbol; logg = 4.0, initial_turbulent_velocity = 2.0,
                                  log_metal = 0.0, file_name = (@__DIR__) * "/phoenix.dat")
    table = readdlm(file_name)
    bit_vector = @. table[:, 1] == initial_turbulent_velocity &&
                    table[:, 3] == logg &&
                    table[:, 5] == log_metal
    table = table[bit_vector, :]
    table = reshape(table, (4, :, 17))

    channels = Dict(
        :u => 6,
        :v => 7,
        :b => 8,
        :y => 9,
        :U => 10,
        :B => 11,
        :V => 12,
        :R => 13,
        :I => 14,
        :J => 15,
        :H => 16,
        :K => 17
    )

    channel_index = channels[channel]

    temperatures :: Vector{Float64} = table[1, :, 4]
    bit_vector = temperatures .<= 6000

    temperatures = temperatures[bit_vector]
    temperatures_range = temperatures[1] : (temperatures[2] - temperatures[1]) : temperatures[end]

    a1 :: Vector{Float64} = table[1, bit_vector, channel_index]
    a2 :: Vector{Float64} = table[2, bit_vector, channel_index]
    a3 :: Vector{Float64} = table[3, bit_vector, channel_index]
    a4 :: Vector{Float64} = table[4, bit_vector, channel_index]

    interpolants = cubic_spline_interpolation.(
        Ref(temperatures_range),
        (a1, a2, a3, a4),
        extrapolation_bc = Flat()
    )

    return x -> x .|> interpolants
end

J_coefs_interpolant = claret_coefs_interpolant(:J)
K_coefs_interpolant = claret_coefs_interpolant(:K)

StructTypes.StructType(::typeof(J_coefs_interpolant)) = StructTypes.StringType()
# StructTypes.StructType(::typeof(K_coefs_interpolant)) = StructTypes.StringType()


end