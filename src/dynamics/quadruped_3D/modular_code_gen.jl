# Structure definition
abstract type AbstractMethods end
abstract type AbstractModules end

mutable struct Quad36Methods <: AbstractMethods
	L::Any
	B::Any
	A::Any
	J::Any
	ϕ::Any
	k::Any
	p_torso::Any
	p_shoulder_1::Any
	p_shoulder_2::Any
	p_shoulder_3::Any
	p_shoulder_4::Any
	p_thigh_1::Any
	p_thigh_2::Any
	p_thigh_3::Any
	p_thigh_4::Any
	p_calf_1::Any
	p_calf_2::Any
	p_calf_3::Any
	p_calf_4::Any
	J_torso::Any
	J_shoulder_1::Any
	J_shoulder_2::Any
	J_shoulder_3::Any
	J_shoulder_4::Any
	J_thigh_1::Any
	J_thigh_2::Any
	J_thigh_3::Any
	J_thigh_4::Any
	J_calf_1::Any
	J_calf_2::Any
	J_calf_3::Any
	J_calf_4::Any
end

function Quad36Methods()
	n_fields = length(Quad36Methods.types)
	mth = Quad36Methods(fill(nothing, n_fields)...)
	return mth
end

struct Quad36Modules{T} <: AbstractModules
	L::AbstractVector{T}
	B::AbstractMatrix{T}
	A::AbstractMatrix{T}
	J::AbstractMatrix{T}
	ϕ::AbstractVector{T}
	k::AbstractVector{T}
	p_torso::AbstractVector{T}
	p_shoulder_1::AbstractVector{T}
	p_shoulder_2::AbstractVector{T}
	p_shoulder_3::AbstractVector{T}
	p_shoulder_4::AbstractVector{T}
	p_thigh_1::AbstractVector{T}
	p_thigh_2::AbstractVector{T}
	p_thigh_3::AbstractVector{T}
	p_thigh_4::AbstractVector{T}
	p_calf_1::AbstractVector{T}
	p_calf_2::AbstractVector{T}
	p_calf_3::AbstractVector{T}
	p_calf_4::AbstractVector{T}
	J_torso::AbstractMatrix{T}
	J_shoulder_1::AbstractMatrix{T}
	J_shoulder_2::AbstractMatrix{T}
	J_shoulder_3::AbstractMatrix{T}
	J_shoulder_4::AbstractMatrix{T}
	J_thigh_1::AbstractMatrix{T}
	J_thigh_2::AbstractMatrix{T}
	J_thigh_3::AbstractMatrix{T}
	J_thigh_4::AbstractMatrix{T}
	J_calf_1::AbstractMatrix{T}
	J_calf_2::AbstractMatrix{T}
	J_calf_3::AbstractMatrix{T}
	J_calf_4::AbstractMatrix{T}
end

function Quad36Modules(dims::Dimensions, env::Environment, T::DataType=Float64)
	nu = dims.u
	nw = dims.w
	nc = dims.c
	np = dim(env)
	L = zeros(T, 1)
	B = zeros(T, nu, nq)
	A = zeros(T, nw, nq)
	J = zeros(T, np * nc, nq)
	ϕ = zeros(T, nc)
	k = zeros(T, np * nc)
	p_torso = zeros(T, np)
	p_shoulder_1 = zeros(T, np)
	p_shoulder_2 = zeros(T, np)
	p_shoulder_3 = zeros(T, np)
	p_shoulder_4 = zeros(T, np)
	p_thigh_1 = zeros(T, np)
	p_thigh_2 = zeros(T, np)
	p_thigh_3 = zeros(T, np)
	p_thigh_4 = zeros(T, np)
	p_calf_1 = zeros(T, np)
	p_calf_2 = zeros(T, np)
	p_calf_3 = zeros(T, np)
	p_calf_4 = zeros(T, np)
	J_torso = zeros(T, np, nq)
	J_shoulder_1 = zeros(T, np, nq)
	J_shoulder_2 = zeros(T, np, nq)
	J_shoulder_3 = zeros(T, np, nq)
	J_shoulder_4 = zeros(T, np, nq)
	J_thigh_1 = zeros(T, np, nq)
	J_thigh_2 = zeros(T, np, nq)
	J_thigh_3 = zeros(T, np, nq)
	J_thigh_4 = zeros(T, np, nq)
	J_calf_1 = zeros(T, np, nq)
	J_calf_2 = zeros(T, np, nq)
	J_calf_3 = zeros(T, np, nq)
	J_calf_4 = zeros(T, np, nq)

	tmp = Quad36Modules{T}(L, B, A, J, ϕ, k,
		p_torso,
		p_shoulder_1,
		p_shoulder_2,
		p_shoulder_3,
		p_shoulder_4,
		p_thigh_1,
		p_thigh_2,
		p_thigh_3,
		p_thigh_4,
		p_calf_1,
		p_calf_2,
		p_calf_3,
		p_calf_4,
		J_torso,
		J_shoulder_1,
		J_shoulder_2,
		J_shoulder_3,
		J_shoulder_4,
		J_thigh_1,
		J_thigh_2,
		J_thigh_3,
		J_thigh_4,
		J_calf_1,
		J_calf_2,
		J_calf_3,
		J_calf_4,
		)
	return tmp
end

struct Quad36{T} <: ContactModel
	dim::Dimensions          # dimensions
	g::T                     # parameters
	μ_world::T
	μ_joint::T
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

	tmp::Quad36Modules{T}    # Temporary variables allowing for non-allocating modular code
	mth::Quad36Methods       # In place (non-allocating) modular methods
	mth_alloc::Quad36Methods # Out of place (allocating) modular methods
end

function generate_modular_expression(model::Quad36, fct::Any, fct_name::Symbol,
	input_size::Vector{Int}, output_size::Vector{Int}; scalar::Bool=false)
	# Declare variables
	num_inputs = length(input_size)
	X = [[Num(Variable(Symbol("X$i"), j)) for j=1:input_size[i]] for i=1:num_inputs]

	# Create, simplify and reshape the output
	y = fct(model, X...)
	y = Symbolics.simplify(y)
	y = reshape(y, output_size...)
	scalar && (y = [y])

	# Build function
	expr = Dict{Symbol, Expr}()
	expr_alloc = Dict{Symbol, Expr}()
	expr_alloc[fct_name], expr[fct_name] = build_function(y, X...)
	return expr, expr_alloc
end

function instantiate_expression(model::ContactModel, expr::Dict{Symbol, Expr},
		expr_alloc::Dict{Symbol, Expr})

	ks = keys(expr)
	@assert ks == keys(expr_alloc)
	for k in ks
		setfield!(model.mth, k,  eval(expr[k]))
		setfield!(model.mth_alloc, k,  eval(expr_alloc[k]))
	end
	return nothing
end




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
g = 9.81      # gravity
μ_world = 1.0 # coefficient of friction
μ_joint = 0.1 # coefficient of torque friction at the joints

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
tmp = Quad36Modules(dims, environment, T)
mth = Quad36Methods()
mth_alloc = Quad36Methods()
quad = Quad36(dims, g, μ_world, μ_joint,
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

expr, expr_alloc = generate_modular_expression(quad, lagrangian, :L, [nq, nq], [1], scalar=true)
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, lagrangian, :L, [nq, nq], [1], scalar=true)
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, B_func, :B, [nq], [nu, nq])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, A_func, :A, [nq], [nw, nq])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, J_func, :J, [nq], [np * nc, nq])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, ϕ_func, :ϕ, [nq], [nc])
instantiate_expression(quad, expr, expr_alloc)
expr, expr_alloc = generate_modular_expression(quad, kinematics, :k, [nq], [np * nc])
instantiate_expression(quad, expr, expr_alloc)


B = zeros(T, nu, nq)
A = zeros(T, nw, nq)
J = zeros(T, np, np * nc)
ϕ = zeros(T, nc)
k = zeros(T, np * nc)
p_torso = zeros(T, np)
p_shoulder_1 = zeros(T, np)
p_shoulder_2 = zeros(T, np)
p_shoulder_3 = zeros(T, np)
p_shoulder_4 = zeros(T, np)
p_thigh_1 = zeros(T, np)
p_thigh_2 = zeros(T, np)
p_thigh_3 = zeros(T, np)
p_thigh_4 = zeros(T, np)
p_calf_1 = zeros(T, np)
p_calf_2 = zeros(T, np)
p_calf_3 = zeros(T, np)
p_calf_4 = zeros(T, np)
J_torso = zeros(T, np, nq)
J_shoulder_1 = zeros(T, np, nq)
J_shoulder_2 = zeros(T, np, nq)
J_shoulder_3 = zeros(T, np, nq)
J_shoulder_4 = zeros(T, np, nq)
J_thigh_1 = zeros(T, np, nq)
J_thigh_2 = zeros(T, np, nq)
J_thigh_3 = zeros(T, np, nq)
J_thigh_4 = zeros(T, np, nq)
J_calf_1 = zeros(T, np, nq)
J_calf_2 = zeros(T, np, nq)
J_calf_3 = zeros(T, np, nq)
J_calf_4 = zeros(T, np, nq)

# Fast - Allocation-free - Modular methods
function L_fast!(model::Quad36{T}, q::AbstractVector{T}, q̇::AbstractVector{T}) where {T}
	model.mth.L(model.tmp.L, q, q̇)
	return nothing
end

function B_fast!(model::Quad36{T}, q::AbstractVector{T}) where {T}
	model.mth.B(model.tmp.B, q)
	return nothing
end

function A_fast!(model::Quad36{T}, q::AbstractVector{T}) where {T}
	model.mth.A(model.tmp.A, q)
	return nothing
end

function J_fast!(model::Quad36{T}, q::AbstractVector{T}) where {T}
	model.mth.J(model.tmp.J, q)
	return nothing
end

function ϕ_fast!(model::Quad36{T}, q::AbstractVector{T}) where {T}
	model.mth.ϕ(model.tmp.ϕ, q)
	return nothing
end

function kinematics_fast!(model::Quad36{T}, q::AbstractVector{T}) where {T}
	model.mth.k(model.tmp.k, q)
	return nothing
end

# Evaluate
q_ = rand(nq)
q̇_ = rand(nq)
L_fast!(quad, q_, q̇_)
@benchmark L_fast!(quad, q_, q̇_)
# @benchmark B_fast!(quad, q_)
# @benchmark A_fast!(quad, q_)
# @benchmark J_fast!(quad, q_)
# @benchmark ϕ_fast!(quad, q_)
@benchmark kinematics_fast!(quad, q_)

@benchmark lagrangian(quad, q_, q̇_)
# @benchmark B_func(quad, q_)
# @benchmark A_func(quad, q_)
# @benchmark J_func(quad, q_)
# @benchmark ϕ_func(quad, q_)
# @benchmark kinematics(quad, q_)


@variables q[1:nq]
@variables q̇[1:nq]

lexpr = lagrangian(quad, q, q̇)
@show lexpr



@variables xx[1:10]
yy = zeros(eltype(xx), 5,5)
yy[1,1] = xx[1]
yy = Symbolics.simplify(reshape(yy, [5,5]...))

@variables q[1:nq]
@variables q̇[1:nq]

# Lagrangian
L = lagrangian(quad, q, q̇)
# L = Symbolics.simplify.(L)
@show L

# dLq = Symbolics.gradient(L, q, simplify=true)
dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
ddL = Symbolics.sparsehessian(L, [q; q̇], simplify=true)
ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

L(q, q̇)
L()
