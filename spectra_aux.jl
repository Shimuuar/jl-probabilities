module spectra_aux

export read_all, correct_columns, read_xy, plot_ribbon


using DelimitedFiles
using Plots
using DataFrames
using Measurements: value, uncertainty

function read_all(filename)
	file = readdlm(filename)
	headers = file[6, :]
	all_data = Float64.(file[7:end, :])
	return DataFrame(all_data, headers)
end

function correct_columns(df_src)
	DataFrame(x = df_src.Uset, y = df_src.CR, y_std = df_src.CRerr)
end

read_xy = correct_columns ∘ read_all


function plot_ribbon(x, y; kwargs...)
    plot(x, value.(y), ribbon = uncertainty.(y), kwargs...)
end

end