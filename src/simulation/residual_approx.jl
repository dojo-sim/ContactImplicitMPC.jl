
function res_con(model::ContactModel, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	[s1 .- ϕ_func(model, env, q2);
	 s2 .- (μ[1] * γ1 .- E_func(model, env) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

function rz_approx!(s, rz, z, θ)
	model = s.model
	env = s.env
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.c * friction_dim(env)

	rz .= 0.0

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)

	# Dynamics
	rz[1:nq, 1:nq] = model.dyn.∂q2(h, q0, q1, u1, w1, λ1, q2) + s.con.dJ(λ1, q2)
	rz[1:nq, 1:(nq + nc + nb)] += s.con.dλ1(q2) * s.con.dcf(γ1, b1, q2, k)

	# Maximum dissipation
	rz[nq + nc .+ (1:nb), nq + nc + nb .+ [1:nc; 2nc .+ (1:nb)]] = s.con.mdψη(vT, ψ1, η1)
	rz[nq + nc .+ (1:nb), 1:nq] += s.con.mdvs(vT, ψ1, η1) * s.con.vsq2(q1, q2, k, h)

	# Other constraints
	s.con.rcz(view(rz, collect([(nq .+ (1:nc))..., ((nq + nc + nb + 1):num_var(model, env))...]), :), z, θ)
end

function rθ_approx!(s, rθ, z, θ)
	model = s.model
	env = s.env
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.c * friction_dim(env)

	rθ .= 0.0

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, env, γ1, b1, q2, k)
	vT = velocity_stack(model, env, q1, q2, k, h)

	# Dynamics
	idx = collect([(1:2nq + model.dim.u + model.dim.w)..., num_data(model)])
	rθ[1:nq, idx] = model.dyn.dθ(h, q0, q1, u1, w1, λ1, q2)

	# Maximum dissipation
	idx = collect([(nq .+ (1:nq))..., num_data(model)])
	rθ[nq + nc .+ (1:nb), idx] = s.con.vsq1h(q1, q2, k, h)

	# Other constraints
	s.con.rcθ(view(rθ, collect([(nq .+ (1:nc))..., ((nq + nc + nb + 1):num_var(model, env))...]), :), z, θ)
end
