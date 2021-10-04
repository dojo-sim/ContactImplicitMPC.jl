
################################################################################
# Model
################################################################################

mutable struct WallPendulum{T} <: ContactModel
	n::Int # state dim
	m::Int # control dim
	mp::T # mass
	l::T # length
	g::T # gravity
	k::T # spring
	d::T # distance to wall
end

function dynamics_model(model::WallPendulum{T}, dt::T; mode::Symbol = :none) where {T}
	mp = model.mp
	l = model.l
	g = model.g
	k = model.k
	d = model.d
	B = dt * [0 1 / (mp*l^2)]'
	if mode == :none
		A = I + dt * [0   1;
			          g/l 0]
		c = dt * [0, 0]
	elseif mode == :left
		A = I + dt * [0          1;
			          g/l - k/mp 0]
		c = dt * [0, k*d / (mp*l)]
	elseif mode == :right
		A = I + dt * [0          1;
			          g/l - k/mp 0]
		c = dt * [0, -k*d / (mp*l)]
	else
		error("Unknown contact mode.")
	end
	return A, B, c
end

function get_mode(x::AbstractVector{T}, u::AbstractVector{T}, model::WallPendulum{T}) where {T}
	d = model.d
	l = model.l
	if -d/l <= x[1] <= d/l # none
		return 1
	elseif x[1] > d/l # left
		return 2
	elseif x[1] < -d/l # right
		return 3
	end
end

function dynamics(x0::AbstractVector{T}, u0::AbstractVector{T},
		A::Vector{Matrix{T}}, B::Vector{Matrix{T}}, c::Vector{Vector{T}},
		model::WallPendulum{T}) where {T}
	i = get_mode(x0, u0, model)
	x1 = A[i] * x0 + B[i] * u0 + c[i]
	return x1
end


function kinematics(x::AbstractVector{T}) where {T}
	θ = x[1] + π/2
	θd = x[2]
	p = l * [cos(θ), sin(θ)]
	return p
end

################################################################################
# Domain
################################################################################

"""
	Box domain for the state-control vector (x, u).
	Domain C = {x_min <= x <= x_max} × {u_min <= u <= u_max}.
	C = {S x + R u <= T}
"""
mutable struct BoxDomain12{T}
	n::Int # state dim
	m::Int # control dim
	x_min::Vector{T} # box boundary
	x_max::Vector{T} # box boundary
	u_min::Vector{T} # box boundary
	u_max::Vector{T} # box boundary
	S::Matrix{T} # affine inequality representation
	R::Matrix{T} # affine inequality representation
	T::Vector{T} # affine inequality representation
end

function BoxDomain12(x_min, x_max, u_min, u_max)
	n = length(x_min)
	m = length(u_min)
	S = Matrix([-I(n); I(n); zeros(m, n); zeros(m, n)])
	R = Matrix([zeros(n, m); zeros(n, m); -I(m); I(m)])
	T = [-x_min; x_max; -u_min; u_max]
	return BoxDomain12(n, m, x_min, x_max, u_min, u_max, S, R, T)
end

function inclusion(C::BoxDomain12, x, u)
	all(C.S * x + C.R * u .<= C.T)
end

function domain(model::WallPendulum{T}; mode::Symbol = :none) where{T}
	m = model.m
	l = model.l
	g = model.g
	k = model.k
	d = model.d

	u_min = [-4.0]
	u_max = [ 4.0]

	if mode == :none
		# @warn "need to revise when adding right wall."
		x_min = [-d/l, -1.5]
		# x_min = [-d/l*2, -1.5]
		x_max = [ d/l,    1.5]
	elseif mode == :left
		x_min = [ d/l,   -1.5]
		x_max = [ d/l*2,  1.5]
	elseif mode == :right
		x_min = [-d/l*2,   -1.5]
		x_max = [-d/l,  1.5]
	else
		error("Unknown contact mode.")
	end
	return BoxDomain12(x_min, x_max, u_min, u_max)
end

# @testset "Domain" begin
# 	# Test
# 	x_min = [-d/l*2, -1.5]
# 	x_max = [ d/l*2,  1.5]
#
# 	x_min_1 = [-d/l*2, -1.5]
# 	x_max_1 = [ d/l,    1.5]
#
# 	x_min_2 = [ d/l,   -1.5]
# 	x_max_2 = [ d/l*2,  1.5]
#
# 	u_min = [-4.0]
# 	u_max = [ 4.0]
#
# 	C0 = BoxDomain12(x_min, x_max, u_min, u_max)
# 	C1 = BoxDomain12(x_min_1, x_max_1, u_min, u_max)
# 	C2 = BoxDomain12(x_min_2, x_max_2, u_min, u_max)
#
# 	xu1 = [[d/l/2,   -1.0], [-2.0]]
# 	xu2 = [[1.5*d/l,  1.0], [ 2.0]]
# 	xu3 = [[2.1d/l,    0.0], [ 0.0]]
#
# 	@test inclusion(C0, xu1...)
# 	@test inclusion(C0, xu2...)
# 	@test !inclusion(C0, xu3...)
# 	@test inclusion(C1, xu1...)
# 	@test !inclusion(C1, xu2...)
# 	@test !inclusion(C1, xu3...)
# 	@test !inclusion(C2, xu1...)
# 	@test inclusion(C2, xu2...)
# 	@test !inclusion(C2, xu3...)
# end

################################################################################
# Problem
################################################################################

mutable struct WallProblem16{F}
	model::WallPendulum{F}
	T::Int
	x0::Vector{F}
	Q::F
	Qf::F
	R::F
	β::F
	A::Vector{Matrix{F}}
	B::Vector{Matrix{F}}
	c::Vector{Vector{F}}
	C::Vector{BoxDomain12{F}}
end

function build_optimizer(prob::WallProblem16{F}) where {F}
	n = prob.model.n
	m = prob.model.m
	nd = length(C)
	T = prob.T
	x0 = prob.x0
	Q = prob.Q
	Qf = prob.Qf
	R = prob.R
	β = prob.β

	Mstar = [β * ones(2n+2m) for i = 1:nd]
	MU = β * ones(n) # upper bound
	ML = - β * ones(n) # lower bound

	optimizer = Model(Gurobi.Optimizer)
	set_silent(optimizer)

	# variables
	@variable(optimizer, x[1:T+1, 1:n])
	@variable(optimizer, u[1:T, 1:m])
	@variable(optimizer, δ[1:nd, 1:T], Bin)
	@variable(optimizer, z[1:nd, 1:T, 1:n])

	# constraints
	# intitial state
	@constraint(optimizer, initial_state, x[1,:] .== x0);
	# 16.a
	@constraint(optimizer, c16a[i in 1:nd, t in 1:T], C[i].S * x[t,:] + C[i].R * u[t,:] - C[i].T .<= Mstar[i] * (1 - δ[i,t]));
	# 16.b
	@constraint(optimizer, c16b[t in 1:T], sum(δ[:, t]) == 1)
	# 18
	@constraint(optimizer, c18[t in 1:T], x[t+1, :] .== sum([z[i,t,:] for i = 1:nd]))
	# 22.a
	@constraint(optimizer, c22a[i in 1:nd, t in 1:T], z[i,t,:] .<= MU * δ[i,t])
	# 22.b
	@constraint(optimizer, c22b[i in 1:nd, t in 1:T], z[i,t,:] .>= ML * δ[i,t])
	# 22.c
	@constraint(optimizer, c22c[i in 1:nd, t in 1:T], z[i,t,:] .<= A[i] * x[t,:] + B[i] * u[t,:] + c[i] - ML * (1 - δ[i,t]))
	# 22.d
	@constraint(optimizer, c22d[i in 1:nd, t in 1:T], z[i,t,:] .>= A[i] * x[t,:] + B[i] * u[t,:] + c[i] - MU * (1 - δ[i,t]))

	# objective
	@objective(optimizer, Min, Q * sum(x[1:T,:].^2) + Qf * sum(x[T+1,:].^2) + R * sum(u.^2))

	return optimizer, x, u, δ, z
end

function get_control(prob::WallProblem16{T}) where {T}
	optimizer, x, u, δ, z = build_optimizer(prob)
	optimize!(optimizer)
	return value.(u)[1,:]
end

function simulate!(prob::WallProblem16{T}, x0::Vector{T}, H::Int;
		w::Vector{Vector{T}} = Vector{Vector{T}}(), iw::Vector{Int} = Vector{Int}()) where {T}
	x = Vector{Vector{T}}([x0])
	u = Vector{Vector{T}}()
	s = Vector{T}()
	for h = 1:H
		println(h, "/", H)
		prob.x0 = deepcopy(x0)
		s0 = @elapsed begin
			u0 = get_control(prob)
		end
		i = findfirst(x -> x == h, iw)
		(i != nothing) && (u0 += w[i])
		x0 = dynamics(x0, u0, prob.A, prob.B, prob.c, prob.model)
		push!(x, deepcopy(x0))
		push!(u, deepcopy(u0))
		push!(s, s0)
	end
	return x, u, s
end




abstract type ControlPolicy{T}
end

mutable struct OpenLoopPolicy16{T} <: ControlPolicy{T}
	u::Vector{Vector{T}}
	N::Int
	N_sample::Int
	count::Int
	dt::T
end

function OpenLoopPolicy16(u::AbstractVector{V}; N_sample::Int = 1, dt::T = 0.01) where {V,T}
	OpenLoopPolicy16(u, length(u), N_sample, 0, dt)
end

mutable struct Simulator12{T}
	p::ControlPolicy{T}
	x::Vector{Vector{T}}
	u::Vector{Vector{T}}
	N_sample::Int
	dt::T
end

function Simulator12(p::ControlPolicy{T}; dt::T = 0.01) where {T}
	x = Vector{Vector{T}}()
	u = Vector{Vector{T}}()
	Simulator12(p, x, u, p.N_sample, p.dt)
end

function simulate!(sim::Simulator12{T}, x0::Vector{T}, H::Int) where {T}
	N_sample = sim.N_sample
	cnt = 0
	x = deepcopy(x0)
	push!(sim.x, deepcopy(x0))
	for t = 1:H
		x, u = step!(p, x, t)
		push!(sim.x, deepcopy(x))
		push!(sim.u, deepcopy(u))
	end
	return nothing
end


function step!(p::OpenLoopPolicy16{T}, x0::AbstractVector{T}, t::Int) where{T}
	@show t
	@show t%p.N_sample
	((t-1) % p.N_sample) == 0 && (p.count += 1)
	u0 = p.u[p.count]
	x1 = dyna(x0, u0)
	return x1, u0
end
