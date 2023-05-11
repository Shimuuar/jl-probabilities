module TuringTurchinRegularization

export Ω, kernel_matrix, cov_matrix, restored_values, restored_spline, make_Ω2_model, apply_model, make_Ω2_model_cached

using SparseArrays
using LinearAlgebra

using Memoize
using Measurements
using BSplineKit
using Turing


Ω(basis, derivative_order) = galerkin_matrix(basis, (Derivative(derivative_order), Derivative(derivative_order)))


"Kᵢⱼ = ∫K(y, xᵢ) bⱼ(y) dy

(yᵢ after convolution) = Kᵢⱼ * cⱼ, 
where cⱼ are the coefficients of the basis functions bⱼ.
"
function kernel_matrix(kernel_function, x_measured, basis)
    return reduce(hcat, (
        galerkin_projection(y -> kernel_function(y, xᵢ), basis) for xᵢ ∈ x_measured
    ))'
end


function cov_matrix(chains::Chains, append_chains = true)
    if append_chains
        return cov(MCMCChains.to_matrix(chains))
    else
        return cov.(MCMCChains.to_vector_of_matrices(chains))
    end
end


function restored_values(chains::Chains, basis::AbstractBSplineBasis{k, T}, eval_at::AbstractArray = collocation_points(basis)) where {k, T}
    coef_samples = Turing.group(chains, :coefs)
    coef_mean = mean(coef_samples)[:, :mean]
    coef_cov = cov_matrix(coef_samples)

    col_mat = collocation_matrix(basis, eval_at, SparseMatrixCSC{T})

    values = col_mat * coef_mean
    values_cov = col_mat * coef_cov * col_mat'
    values_std = sqrt.(diag(values_cov))

    return values .± values_std
end

function restored_values(chains::Chains, model::Turing.Model)
    basis = model.args.spline_basis
    restored_values(chains, basis)
end
function restored_values(chains, model::Turing.Model, eval_at::AbstractArray)
    basis = model.args.spline_basis
    restored_values(chains, basis, eval_at)
end


function restored_spline(chains::Chains, basis::AbstractBSplineBasis)
    coef_samples = Turing.group(chains, :coefs)
    coef_mean = mean(coef_samples)[:, :mean]
    return Spline(basis, coef_mean)
end


function make_Ω2_model(df, model_func, kernel_func; spline_knots)
    b = BSplineBasis(4, copy(spline_knots))
    b = RecombinedBSplineBasis(Derivative(0), b)

    K = kernel_matrix(kernel_func, df.x, b)
    Ω2 = Ω(b, 2)

    return model_func(df; Ω = Ω2, kernel_matrix = K, kernel_func = kernel_func, spline_basis = b)
end

@memoize Dict function make_Ω2_model_cached(df, model_expr, kernel_func; spline_knots)
    b = BSplineBasis(4, copy(spline_knots))
    b = RecombinedBSplineBasis(Derivative(0), b)

    K = kernel_matrix(kernel_func, df.x, b)
    Ω2 = Ω(b, 2)

    model_func = eval(:(@model $model_expr))

    return model_func(df; Ω = Ω2, kernel_matrix = K, kernel_func = kernel_func, spline_basis = b)
end


function apply_model(model::Turing.Model, basis::AbstractBSplineBasis, sampler)
    samples = sample(model, sample_kwargs)
end



# Convinience functions that are defined for Chains but not ChainDataFrame
# PR to MCMCChains.jl?

"""
    namesingroup(chain_df::::ChainDataFrame, sym::Union{AbstractString,Symbol}; index_type::Symbol=:bracket)
Similarly to `namesingroup(chains::Chains, ...)`,
return the parameters with the same name `sym`, but have a different index. Bracket indexing format
in the form of `:sym[index]` is assumed by default. Use `index_type=:dot` for parameters with dot 
indexing, i.e. `:sym.index`.
If the chain contains a parameter of name `:sym` it will be returned as well.
# Example
```jldoctest
julia> chain_df = ChainDataFrame((
	parameters = Symbol.(["A[1]", "A[2]"]),
	values = [1, 2]
));
julia> namesingroup(chain_df, :A)
2-element Vector{Symbol}:
 Symbol("A[1]")
 Symbol("A[2]")
```
"""
Turing.namesingroup(chain_df::ChainDataFrame, sym::AbstractString; kwargs...) = namesingroup(chain_df, Symbol(sym); kwargs...)
function namesingroup(chain_df::ChainDataFrame, sym::Symbol; index_type::Symbol=:bracket)
    if index_type !== :bracket && index_type !== :dot
        error("index_type must be :bracket or :dot")
    end
    idx_str = index_type == :bracket ? "[" : "."
    # Start by looking up the symbols in the list of parameter names.
    names_of_params = chain_df[:, :parameters]
    regex = Regex("^\\Q$sym\\E\$|^\\Q$sym$idx_str\\E")
    indices = findall(x -> match(regex, string(x)) !== nothing, names_of_params)
    return names_of_params[indices]
end

"""
    group(chain_df::ChainDataFrame, name::Union{AbstractString,Symbol}; index_type::Symbol=:bracket)
Similarly to `group(chains::Chains, ...)`,
Return a subset of the ChainDataFrame containing parameters with the same `name`, but a different index.
Bracket indexing format in the form of `:name[index]` is assumed by default. Use `index_type=:dot` for parameters with dot 
indexing, i.e. `:sym.index`.
"""
function Turing.group(chain_df::ChainDataFrame, name::Union{AbstractString,Symbol}; kwargs...)
    return chain_df[namesingroup(chain_df, name; kwargs...)]
end

end
