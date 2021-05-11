# """
# 	generate_base_expressions(model::ContactDynamicsModel)
# Generate fast base methods using Symbolics symbolic computing tools.
# """
# function generate_modular_base_expressions(model::ContactDynamicsModel;
# 	M_analytical = true)
#
# 	nq = model.dim.q
# 	nu = model.dim.u
# 	nw = model.dim.w
# 	nc = model.dim.c
# 	nb = model.dim.b
# 	np = dim(model.env)
#
# 	# Declare variables
# 	@variables q[1:nq]
# 	@variables q̇[1:nq]
#
# 	# Lagrangian
# 	L = lagrangian(model, q, q̇)
# 	L = Symbolics.simplify.(L)
#
# 	dLq = Symbolics.gradient(L, q, simplify=true)
# 	dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
# 	ddL = Symbolics.sparsehessian(L, [q; q̇], simplify=true)
# 	ddLq̇q = ddL[nq .+ (1:nq), 1:nq]
#
# 	# Mass Matrix
# 	if M_analytical
# 		M = M_func(model, q)
# 		M = reshape(M, (nq, nq))
# 		M = Symbolics.simplify.(M)
# 	else
# 		M = ddL[nq .+ (1:nq), nq .+ (1:nq)]
# 	end
#
# 	# Control input Jacobian
# 	B = B_func(model, q)
# 	B = reshape(B, (nu, nq))
# 	B = Symbolics.simplify.(B)
#
# 	# Disturbance input Jacobian
# 	A = A_func(model, q)
# 	A = reshape(A, (nw, nq))
# 	A = Symbolics.simplify.(A)
#
# 	# Contact Jacobian
# 	J = J_func(model, q)
# 	J = reshape(J, (np * nc, nq))
# 	J = Symbolics.simplify.(J)
#
# 	# Coriolis and Centrifugal forces Jacobians
# 	C = ddLq̇q * q̇ - dLq
# 	C = Symbolics.simplify.(C)
#
# 	# Kinematics
# 	k = kinematics(model, q)
# 	k = simplify.(k)
#
# 	# Build function
# 	expr = Dict{Symbol, Expr}()
# 	expr[:L]    = build_function([L], q, q̇)[1] # need to transpose to get a line vector
# 	expr[:M]    = build_function(M, q)[1]
# 	expr[:B]    = build_function(B, q)[1]
# 	expr[:A]    = build_function(A, q)[1]
# 	expr[:J]    = build_function(J, q)[1]
# 	expr[:C]    = build_function(C, q, q̇)[1]
# 	expr[:k]    = build_function(k, q)[1]
#
# 	return expr
# end


# function generate_modular_ϕ(model::ContactDynamicsModel)
# 	nq = model.dim.q
# 	nu = model.dim.u
# 	nw = model.dim.w
# 	nc = model.dim.c
# 	nb = model.dim.b
# 	np = dim(model.env)
#
# 	# Declare variables
# 	@variables q[1:nq]
#
# 	# Lagrangian
# 	ϕ = Symbolics.simplify(ϕ_func(model, q))
#
#
# 	# Build function
# 	expr = Dict{Symbol, Expr}()
# 	expr[:L]    = build_function([L], q, q̇)[1] # need to transpose to get a line vector
# 	expr[:M]    = build_function(M, q)[1]
# 	return expr
# end


# model = quadruped_3D
# nq = model.dim.q
# nu = model.dim.u
# nw = model.dim.w
# nc = model.dim.c
# nb = model.dim.b
# np = dim(model.env)
#
# # Declare variables
# @variables q[1:nq]
#
# ϕ = ϕ_func(model, q)
# ϕ = Symbolics.simplify(ϕ)
# ϕfast1 = eval(build_function(ϕ, q)[1])
# ϕfast2 = eval(build_function(ϕ, q)[2])
#
# q_ = rand(nq)
# ϕ1_ = rand(nc)
# ϕ2_ = rand(nc)
# @benchmark ϕ1_ = ϕfast1(q_)
# @benchmark ϕfast2(ϕ2_, q_)


# function generate_modular_expressions(model::Quad32)
# 	nq = model.dim.q
# 	nu = model.dim.u
# 	nw = model.dim.w
# 	nc = model.dim.c
# 	nb = model.dim.b
#
# 	# Declare variables
# 	@variables q[1:nq]
#
# 	# Collision
# 	ϕ = ϕ_func(model, q)
# 	ϕ = Symbolics.simplify(ϕ)
#
# 	# Build function
# 	expr = Dict{Symbol, Expr}()
# 	expr_alloc = Dict{Symbol, Expr}()
#
# 	expr[:ϕ] = build_function(ϕ, q)[1]
# 	expr_alloc[:ϕ] = build_function(ϕ, q)[2]
# 	return expr, expr_alloc
# end


# Structure definition
abstract type AbstractMethods end
abstract type AbstractModules end

mutable struct Quad32Methods <: AbstractMethods
	B::Any
	A::Any
	J::Any
	ϕ::Any
	k::Any
end

function Quad32Methods()
	n_fields = length(Quad32Methods.types)
	mth = Quad32Methods(fill(nothing, n_fields)...)
	return mth
end

struct Quad32Modules{T} <: AbstractModules
	B::AbstractMatrix{T}
	A::AbstractMatrix{T}
	J::AbstractMatrix{T}
	ϕ::AbstractVector{T}
	k::AbstractVector{T}
end

function Quad32Modules(dims::Dimensions, env::Environment, T::DataType=Float64)
	nu = dims.u
	nw = dims.w
	nc = dims.c
	np = dim(env)
	B = zeros(T, nu, nq)
	A = zeros(T, nw, nq)
	J = zeros(T, np * nc, nq)
	ϕ = zeros(T, nc)
	k = zeros(T, np * nc)
	tmp = Quad32Modules{T}(B, A, J, ϕ, k, )
	return tmp
end

struct Quad32{T} <: ContactDynamicsModel
	dim::Dimensions          # dimensions
	μ::T                     # parameters
	orientation::Symbol
	shoulder_lateral_offset::T

	# torso
	l_torso::T
	d_torso::T
	m_torso::T
	J_torso::T

	# leg 1
		# shoulder
	l_shoulder1::T
	d_shoulder1::T
	m_shoulder1::T
	J_shoulder1::T
		# thigh
	l_thigh1::T
	d_thigh1::T
	m_thigh1::T
	J_thigh1::T
		# calf
	l_calf1::T
	d_calf1::T
	m_calf1::T
	J_calf1::T

	# leg 2
		# shoulder
	l_shoulder2::T
	d_shoulder2::T
	m_shoulder2::T
	J_shoulder2::T
		# thigh
	l_thigh2::T
	d_thigh2::T
	m_thigh2::T
	J_thigh2::T
		# calf
	l_calf2::T
	d_calf2::T
	m_calf2::T
	J_calf2::T

	# leg 3
		# shoulder
	l_shoulder3::T
	d_shoulder3::T
	m_shoulder3::T
	J_shoulder3::T
		# thigh
	l_thigh3::T
	d_thigh3::T
	m_thigh3::T
	J_thigh3::T
		# calf
	l_calf3::T
	d_calf3::T
	m_calf3::T
	J_calf3::T

	# leg 4
		# shoulder
	l_shoulder4::T
	d_shoulder4::T
	m_shoulder4::T
	J_shoulder4::T
		# thigh
	l_thigh4::T
	d_thigh4::T
	m_thigh4::T
	J_thigh4::T
		# calf
	l_calf4::T
	d_calf4::T
	m_calf4::T
	J_calf4::T

	env::Environment

	tmp::Quad32Modules{T}    # Temporary variables allowing for non-allocating modular code
	mth::Quad32Methods       # In place (non-allocating) modular methods
	mth_alloc::Quad32Methods # Out of place (allocating) modular methods
end

function generate_modular_expression(model::Quad32, fct::Any, fct_name::Symbol,
	input_size::Int, output_size::Vector{Int})
	# Declare variables
	@variables x[1:input_size]

	# Function
	y = fct(model, x)
	y = Symbolics.simplify(reshape(y, output_size...))

	# Build function
	expr = Dict{Symbol, Expr}()
	expr_alloc = Dict{Symbol, Expr}()

	expr_alloc[fct_name], expr[fct_name] = build_function(y, x)
	return expr, expr_alloc
end

function instantiate_expression(model::ContactDynamicsModel, expr::Dict{Symbol, Expr},
		expr_alloc::Dict{Symbol, Expr})

	ks = keys(expr)
	@assert ks == keys(expr_alloc)
	for k in ks
		setfield!(model.mth, k,  eval(expr[k]))
		setfield!(model.mth_alloc, k,  eval(expr_alloc[k]))
	end
	return nothing
end

@variables xx[1:10]
yy = zeros(eltype(xx), 5,5)
yy[1,1] = xx[1]
yy = Symbolics.simplify(reshape(yy, [5,5]...))



# Dimensions
environment = environment_3D_flat()
nq = 18
nu = 12
nw = 3
nc = 4
nb = nc*4
np = dim(environment)
dims = Dimensions(nq, nu, nw, nc, nb)

# Parameters
μ_ = 1.0
shoulder_lateral_offset = 0.047
orientation = :MRP
m_torso = 4.713
m_shoulder = 0.696
m_thigh = 1.013
m_leg = 0.166

J_torso = 0.01683
J_shoulder = 0.000469246
J_thigh = 0.00552
J_leg = 0.00299

l_torso = 0.183*2
l_shoulder = 0.08505
l_thigh = 0.2
l_leg = 0.2

d_torso = 0.5 * l_torso + 0.0127
d_shoulder = 0.5 * l_shoulder -0.003311
d_thigh = 0.5 * l_thigh - 0.00323
d_leg = 0.5 * l_leg - 0.006435


# Helpers
T=  Float64
tmp = Quad32Modules(dims, environment, T)
mth = Quad32Methods()
mth_alloc = Quad32Methods()
quad = Quad32(dims, μ_,
	orientation, shoulder_lateral_offset,
	l_torso, d_torso, m_torso, J_torso,
	l_shoulder, d_shoulder, m_shoulder, J_shoulder,
	l_thigh, d_thigh, m_thigh, J_thigh,
	l_leg, d_leg, m_leg, J_leg,
	l_shoulder, d_shoulder, m_shoulder, J_shoulder,
	l_thigh, d_thigh, m_thigh, J_thigh,
	l_leg, d_leg, m_leg, J_leg,
	l_shoulder, d_shoulder, m_shoulder, J_shoulder,
	l_thigh, d_thigh, m_thigh, J_thigh,
	l_leg, d_leg, m_leg, J_leg,
	l_shoulder, d_shoulder, m_shoulder, J_shoulder,
	l_thigh, d_thigh, m_thigh, J_thigh,
	l_leg, d_leg, m_leg, J_leg,
	environment,
	tmp, mth, mth_alloc)

include("modular_model.jl")

expr, expr_alloc = generate_modular_expression(quad, B_func, :B, nq, [nu, nq])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, A_func, :A, nq, [nw, nq])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, J_func, :J, nq, [np * nc, nq])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, ϕ_func, :ϕ, nq, [nc])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, kinematics, :k, nq, [np * nc])
instantiate_expression(quad, expr, expr_alloc)


B = zeros(T, nu, nq)
A = zeros(T, nw, nq)
J = zeros(T, np, np * nc)
ϕ = zeros(T, nc)
k = zeros(T, np * nc)


# Fast - Allocation-free - Modular methods
function B_fast!(model::Quad32{T}, q::AbstractVector{T}) where {T}
	model.mth.B(model.tmp.B, q)
	return nothing
end

function A_fast!(model::Quad32{T}, q::AbstractVector{T}) where {T}
	model.mth.A(model.tmp.A, q)
	return nothing
end

function J_fast!(model::Quad32{T}, q::AbstractVector{T}) where {T}
	model.mth.J(model.tmp.J, q)
	return nothing
end

function ϕ_fast!(model::Quad32{T}, q::AbstractVector{T}) where {T}
	model.mth.ϕ(model.tmp.ϕ, q)
	return nothing
end

function kinematics_fast!(model::Quad32{T}, q::AbstractVector{T}) where {T}
	model.mth.k(model.tmp.k, q)
	return nothing
end

# Evaluate
q_ = rand(nq)
@benchmark B_fast!(quad, q_)
@benchmark A_fast!(quad, q_)
@benchmark J_fast!(quad, q_)
@benchmark ϕ_fast!(quad, q_)
@benchmark kinematics_fast!(quad, q_)

@benchmark B_func(quad, q_)
@benchmark A_func(quad, q_)
@benchmark J_func(quad, q_)
@benchmark ϕ_func(quad, q_)
@benchmark kinematics(quad, q_)

J_func(quad, q_)
