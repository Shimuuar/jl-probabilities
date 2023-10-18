module ChainCache

using Serialization
using SHA

using StructTypes
using JSON3
using Turing

include("LuminocityModels.jl")
using .LuminocityModels

export
    cached_sample

StructTypes.StructType(::Type{MeshParams}) = StructTypes.Struct()
StructTypes.StructType(::Type{ModelParams}) = StructTypes.Struct()
StructTypes.StructType(::Type{ChainParams}) = StructTypes.Struct()

function cached_sample(chain_params, cache_dir::String="cache")
    chain_params_hash = chain_params |> JSON3.write |> sha1 |> bytes2hex

    cache_dir = joinpath(cache_dir, chain_params_hash)
    cache_file = joinpath(cache_dir, "chain.jls")
    if isfile(cache_file)
        return deserialize(cache_file)
    else
        mkpath(cache_dir)
        samples = _uncached_sample(chain_params)
        serialize(cache_file, samples)

        description_file = joinpath(cache_dir, "parameters.json")
        open(description_file, "w") do f
            write(f, good_json(chain_params))
        end
        return samples
    end
end


function _uncached_sample(chain_params)
    model = model_from_params(chain_params.model_params)
    sampler = eval(Meta.parse(chain_params.sampler_str))

    return sample(
        model,
        sampler,
        chain_params.n_samples,
        init_params=chain_params.init_params
    )
end



"""
Outputs something like:
{
    "a": 1,
    "b": [ 1, 2, 3 ]
}
"""
function good_json(object)
    buffer = IOBuffer()
    JSON3.pretty(buffer, object)
    s = String(take!(buffer))

    arrays = findall(r"\[[^\]]+\]", s)

    for a in arrays[end:-1:begin]
        compact_array = replace(s[a], r"\s+" => " ")
        s = replace(s, s[a] => compact_array)
    end
    return s
end

end