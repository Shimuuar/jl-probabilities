module ChainCache

using Serialization
using SHA

using JSON3
using Turing

export
    cached_sample,
    custom_pretty_json

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
            write(f, custom_pretty_json(chain_params))
        end
        return samples
    end
end


function _uncached_sample(chain_params)
    model_params = chain_params.model_params
    model = model_params.model_function(model_params)

    if chain_params.init_params !== nothing
        syms = DynamicPPL.syms(DynamicPPL.VarInfo(model))
        init_params_array = chain_params.init_params[syms]
    else
        init_params_array = nothing
    end

    return sample(
        model,
        chain_params.sampler,
        chain_params.n_samples,
        init_params=init_params_array
    )
end



"""
Outputs something like:
{
    "a": 1,
    "b": [1, 2, 3]
}
"""
function custom_pretty_json(object)
    object = JSON3.read(JSON3.write(object))

    function format_item(item, level)
        ind = " " ^ (4 * level)
        if isa(item, AbstractDict)
            return "{\n" * join(["$(ind)    \"$(k)\": $(format_item(v, level + 1))" for (k, v) in item], ",\n") * "\n$(ind)}"
        elseif isa(item, AbstractArray)
            return "[" * join([format_item(v, level + 1) for v in item], ", ") * "]"
        else
            return JSON3.write(item)
        end
    end

    return format_item(object, 0)
end

end