mutable struct Simulation{T,W,FC}
	model::Model{T}
	env::Environment{W,FC}
	con::ContactMethods
	res::ResidualMethods
	rz::Any
	rθ::Any 
	ϕ::Any
end

function Simulation(model::Model{T}, env::Environment) where T
	Simulation(model, env, ContactMethods(), ResidualMethods(), zeros(0, 0), zeros(0, 0), q -> nothing)
end

function get_simulation(model::String, env::String, sim_name::String;
		model_variable_name::String = model,
		dynamics_name::String = "dynamics",
		gen_base = true,
		gen_dyn = true,
		approx = false)

	#TODO: assert model exists

	dir_model = joinpath(@__DIR__, "..", "dynamics", model)
	dir_sim   = joinpath(@__DIR__, "..", "simulation", model, sim_name)

	dir_base = joinpath(dir_model, dynamics_name, "base.jld2")
	dir_dyn = joinpath(dir_model, dynamics_name, "dynamics.jld2")
	dir_res = joinpath(dir_sim, "residual.jld2")
	dir_jac = joinpath(dir_sim, "jacobians.jld2")

	model = deepcopy(eval(Symbol(model_variable_name)))
	env = deepcopy(eval(Symbol(env)))
	sim = Simulation(model, env)

	instantiate_base!(sim.model, dir_base)

	instantiate_dynamics!(sim.model, dir_dyn,
		derivs = approx)

	instantiate_residual!(sim,
		dir_res, dir_jac,
		jacobians = (approx ? :approx : :full))

	# signed-distance function 
	@variables q[1:sim.model.nq] 
	ϕ = ϕ_func(sim.model, sim.env, q) 
	sim.ϕ = eval(Symbolics.build_function(ϕ, q, checkbounds=true)[2])

	return sim
end

function z_initialize!(z, model::Model, env::Environment{<:World,LinearizedCone}, q)
	z .= 1.0
	iq2 = index_q2(model, env, quat = false)
	z[iq2] = q
end

function z_initialize!(z::Vector{T}, idx::SVector{nq2,Int}, q::Vector{T}) where {T,nq2} 
	z .= 1.0
	z[idx] .= q
	return nothing
end

function z_initialize!(z, model::Model, env::Environment{<:World,NonlinearCone}, q1) #TODO: make allocation free
	iq2 = index_q2(model, env, quat = false)
	ib1 = index_b1(model, env, quat = false)
	iψ1 = index_ψ1(model, env, quat = false)
	iη1 = index_η1(model, env, quat = false)
	is2 = index_s2(model, env, quat = false)

    z .= 0.1
    z[iq2] = q1

	# second-order cone initializations
	z[ib1] .= 0.1 # primal cones: vector part # TODO redundant
	z[iψ1] .= 1.0 # primal cones: scalar part
	z[iη1] .= 0.1 # dual cones: vector part # TODO redundant
	z[is2] .= 1.0 # dual cones: scalar part
	return z
end

function z_initialize!(z, idx::SVector{nq2,Int}, model::Model, env::Environment{<:World,NonlinearCone}, q1) where {T,nq2}#TODO: make allocation free
	@warn "soc z initialization not implemented"
end

function z_warmstart!(z, model::Model, env::Environment{<:World,LinearizedCone}, q, a)
	iq2 = index_q2(model, env, quat = false)
	iort = index_ort(model, env, quat = false)
	isoc = index_soc(model, env, quat = false)
	z[iq2] = q
	for i in eachindex(iort)
		z[iort[i]] .+= a * rand(length(iort[i]))
	end
	for i in eachindex(isoc)
		for j in eachindex(isoc[i])
			z[isoc[i][j]] .+= a * rand(length(isoc[i][j]))
		end
	end
	nothing
end

function z_warmstart!(z, model::Model, env::Environment{<:World,NonlinearCone}, q, a)
	# @warn "warm start for second-order cone not implemented"
	z_initialize!(z, model, env, q)
	nothing
end

function θ_initialize!(θ, model::Model, q0, q1, u, w, μ, h)
	nq = model.nq
	nu = model.nu
	nw = model.nw
	off = 0
    θ[1:nq] = q0
	off += nq
    θ[off .+ (1:nq)] = q1
	off += nq
    θ[off .+ (1:nu)] = u
	off += nu
    θ[off .+ (1:nw)] = w
	off += nw
	θ[off .+ (1:1)] .= μ
	off += 1
    θ[off .+ (1:1)] .= h
end

function E_func(model::Model, env::Environment{<:World,LinearizedCone})
	nc = model.nc
	nb = nc * friction_dim(env)
	SMatrix{nc, nb}(kron(Diagonal(ones(nc)), ones(1, Int(nb / nc))))
end

function residual(model::Model, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	nc = model.nc
	nb = nc * friction_dim(env)
	nf = Int(nb / nc)
	np = dim(env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	ϕ = ϕ_func(model, env, q2)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	Λ1 = transpose(J_func(model, env, q2)) * λ1 #@@@@ maybe need to use J_fast
	vT_stack = velocity_stack(model, env, q1, q2, k, h)
	ψ_stack = transpose(E_func(model, env)) * ψ1

	# @warn "define residual order"
	[model.dyn.d(h, q0, q1, u1, w1, Λ1, q2);
	 s1 - ϕ;
	 η1 - vT_stack - ψ_stack;
	 s2 .- (μ[1] * γ1 .- E_func(model, env) * b1);
	 γ1 .* s1 .- κ[1];
	 b1 .* η1 .- κ[1];
	 ψ1 .* s2 .- κ[1]]
end

function residual(model::Model, env::Environment{<:World,NonlinearCone}, z, θ, κ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)
	ϕ = ϕ_func(model, env, q2)
	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)
	nc = model.nc
	nf = friction_dim(env)
	ne = dim(env)

	# @warn "define residual order"
	[dynamics(model, env, h, q0, q1, u1, w1, λ1, q2);
	 s1 - ϕ;
	 vcat([η1[(i - 1) * nf .+ (1:nf)] - vT[(i - 1) * nf .+ (1:nf)] for i = 1:nc]...);
	 s2 - μ[1] * γ1;
	 γ1 .* s1 .- κ[1];
	 vcat([
	 	second_order_cone_product( # [ψ1; η1] ∘ [s2; b1] - [κ; 0...0]
			# [ψ1[i]; b1[(i-1) * nf .+ (1:nf)]],
			[ψ1[i]; η1[(i-1) * nf .+ (1:nf)]],
			# [s2[i]; η1[(i-1) * nf .+ (1:nf)]]
			[s2[i]; b1[(i-1) * nf .+ (1:nf)]]
			) - [κ[1]; zeros(nf)]
		for i = 1:nc]...)
	 ]
end

function ResidualMethods()
	function f()
		error("Not Implemented: use instantiate_residual!")
		return nothing
	end
	return ResidualMethods(fill(f, 3)...)
end
