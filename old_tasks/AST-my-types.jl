### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 33c8313e-12ec-450c-b271-594bc8f3d2d2
abstract type Expression end

# ╔═╡ 427c97b6-b036-4019-98a0-fb7eb2cac1af
struct Lit <:Expression
	val
end

# ╔═╡ bf508c04-1a58-47d4-9900-bf53126fca04
struct Var <: Expression
	name::String
end

# ╔═╡ 9cc0450b-e38a-4e6e-bb95-17f3d81e0785
Base.convert(::Type{Expression}, name::String) = Var(name)

# ╔═╡ 659d9f13-6372-4a2a-8fdd-e94ec46540cc
Base.convert(::Type{Expression}, val::Number) = Lit(val)

# ╔═╡ 01739ad1-2422-4ebc-bdfc-ce9ed8d675a2
struct Operation <: Expression
	head::Symbol
	subexpr1::Expression
	subexpr2::Expression
end

# ╔═╡ 72e3f117-e6ff-439d-b997-fb6ecd440b43
Operation(:*, "x", 12)

# ╔═╡ d912cc35-c1cc-42f3-b082-2e0242bdf2e9
md"### Вычислитель"

# ╔═╡ 7669e12e-e5f1-4562-ad6e-3eea5c6c6593
function substitute_vars(var::Var, dict)
	try
		return Lit(dict[var.name])
	catch e
		if isa(e, KeyError)
			return var
		end
	end
end

# ╔═╡ 23202800-340d-4773-a049-1317c22febc2
substitute_vars(lit::Lit, dict) = lit

# ╔═╡ a678dd3f-3784-475a-baba-58da0412c8e2
substitute_vars(operation::Operation, dict) = Operation(
	operation.head,
	substitute_vars(operation.subexpr1, dict),
	substitute_vars(operation.subexpr2, dict)
)

# ╔═╡ 7e5c26b1-0f34-4878-913e-fa8242795cb0
function fold_constants_nonrecursively(head::Symbol, lit1::Lit, lit2::Lit)
	if head == :+
		return Lit(lit1.val + lit2.val)
	elseif head == :*
		return Lit(lit1.val * lit2.val)
	else
		error("operation $head is not supported")
	end
end

# ╔═╡ 088df5fc-e91f-4df7-b99f-fe9bec455105
fold_constants_nonrecursively(head::Symbol, subexpr1::Expression, subexpr2::Expression) = Operation(head, subexpr1, subexpr2)

# ╔═╡ f9683576-54ca-4b27-a3ed-d30c177133a4
function fold_constants(operation::Operation)
	head = operation.head
	subexpr1 = fold_constants(operation.subexpr1)
	subexpr2 = fold_constants(operation.subexpr2)
	return fold_constants_nonrecursively(head, subexpr1, subexpr2)
end

# ╔═╡ 7bf2e82e-b8cb-4cdd-aa0b-33c358ef66be
fold_constants(lit::Lit) = lit

# ╔═╡ 811a4ac6-bf69-4e69-b581-1a482d9ca038
fold_constants(var::Var) = var

# ╔═╡ 351eae00-5180-441e-afa3-d9c628a03d90
fold_constants(Operation(:+, 1, 2))

# ╔═╡ 45bda75a-7212-4ecb-9df7-f0ba52991e88
evaluate(operation::Operation, dict) = fold_constants(
	substitute_vars(operation, dict)
)

# ╔═╡ 30db42b3-5112-4566-83bd-97446786ccad
md"### Преобразование в строку"

# ╔═╡ a2e1855f-1f15-42ef-81aa-69ef602714a1
function Base.string(operation::Operation)
	subexpr1 = operation.subexpr1
	subexpr2 = operation.subexpr2
	s1 = string(subexpr1)
	s2 = string(subexpr2)
	if operation.head == :*
		if isa(subexpr1, Operation) && subexpr1.head == :+
			s1 = "(" * s1 * ")"
		elseif isa(subexpr2, Operation) && subexpr2.head == :+
			s2 = "(" * s2 * ")"
		end
	end
			
	s1 *  string(operation.head) * s2
end

# ╔═╡ 595cc0f0-4f7f-4c5c-ac2b-223991ae19b0
Base.string(lit::Lit) = string(lit.val)

# ╔═╡ 801e4d65-c7eb-4823-95b1-50e2bbe4e8d7
Base.string(var::Var) = var.name

# ╔═╡ 9f93c27d-ff95-460b-a01e-712c98c25ff7
md"### Тесты"

# ╔═╡ 4fdba3d7-19bc-4fb2-af96-3c5f80f91158
operation1 = Operation(
	:+,
	Operation(:*, "x", 12),
	4
)

# ╔═╡ 59cc8a38-f494-43e5-aaef-f40e9948423f
string(operation1)

# ╔═╡ bd06a452-943b-4b77-b4a4-847e38413ade
dict = Dict("x" => 2)

# ╔═╡ 3b151c3e-8422-47da-a697-5903cbaa2103
substitute_vars(Operation(:+, 1, "x"), dict)

# ╔═╡ db153a4b-89ab-4e48-acb7-73b217473952
evaluate(operation1, dict)

# ╔═╡ 555ff9fa-bd27-48e0-add5-ba0f2556ee48
2*12 + 4

# ╔═╡ 003b5064-fa04-4453-99e1-5343a29e07f7
operation2 = Operation(
	:*,
	Operation(:+, "x", 12),
	Operation(:*, 6, 5),
)

# ╔═╡ 3f3a39d8-eea6-4c0c-802e-773443e541cc
string(operation2)

# ╔═╡ 9f275f44-c0d1-47b5-ab54-983dbc1ddd0c
operation2 |> fold_constants |> string

# ╔═╡ 0b2a6ea5-afed-4f86-8419-15e48208b96e
evaluate(operation2, dict)

# ╔═╡ 16a5ff60-7265-4f57-b57b-f5bef56c63b4
(2+12) * 6 * 5

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═33c8313e-12ec-450c-b271-594bc8f3d2d2
# ╠═427c97b6-b036-4019-98a0-fb7eb2cac1af
# ╠═bf508c04-1a58-47d4-9900-bf53126fca04
# ╠═9cc0450b-e38a-4e6e-bb95-17f3d81e0785
# ╠═659d9f13-6372-4a2a-8fdd-e94ec46540cc
# ╠═01739ad1-2422-4ebc-bdfc-ce9ed8d675a2
# ╠═72e3f117-e6ff-439d-b997-fb6ecd440b43
# ╟─d912cc35-c1cc-42f3-b082-2e0242bdf2e9
# ╠═7669e12e-e5f1-4562-ad6e-3eea5c6c6593
# ╠═23202800-340d-4773-a049-1317c22febc2
# ╠═a678dd3f-3784-475a-baba-58da0412c8e2
# ╠═3b151c3e-8422-47da-a697-5903cbaa2103
# ╠═f9683576-54ca-4b27-a3ed-d30c177133a4
# ╠═7e5c26b1-0f34-4878-913e-fa8242795cb0
# ╠═088df5fc-e91f-4df7-b99f-fe9bec455105
# ╠═7bf2e82e-b8cb-4cdd-aa0b-33c358ef66be
# ╠═811a4ac6-bf69-4e69-b581-1a482d9ca038
# ╠═351eae00-5180-441e-afa3-d9c628a03d90
# ╠═45bda75a-7212-4ecb-9df7-f0ba52991e88
# ╟─30db42b3-5112-4566-83bd-97446786ccad
# ╠═a2e1855f-1f15-42ef-81aa-69ef602714a1
# ╠═595cc0f0-4f7f-4c5c-ac2b-223991ae19b0
# ╠═801e4d65-c7eb-4823-95b1-50e2bbe4e8d7
# ╟─9f93c27d-ff95-460b-a01e-712c98c25ff7
# ╠═4fdba3d7-19bc-4fb2-af96-3c5f80f91158
# ╠═59cc8a38-f494-43e5-aaef-f40e9948423f
# ╠═bd06a452-943b-4b77-b4a4-847e38413ade
# ╠═db153a4b-89ab-4e48-acb7-73b217473952
# ╠═555ff9fa-bd27-48e0-add5-ba0f2556ee48
# ╠═003b5064-fa04-4453-99e1-5343a29e07f7
# ╠═3f3a39d8-eea6-4c0c-802e-773443e541cc
# ╠═9f275f44-c0d1-47b5-ab54-983dbc1ddd0c
# ╠═0b2a6ea5-afed-4f86-8419-15e48208b96e
# ╠═16a5ff60-7265-4f57-b57b-f5bef56c63b4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
