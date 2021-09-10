
function res_con(model::ContactModel, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	@warn "define residual order"
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
	idyn = index_dyn(model, env, quat = false)
	iimp = index_imp(model, env, quat = false)
	imdp = index_mdp(model, env, quat = false)
	ifri = index_fri(model, env, quat = false)
	ibimp = index_bimp(model, env, quat = false)
	ibmdp = index_bmdp(model, env, quat = false)
	ibfri = index_bfri(model, env, quat = false)

	iq2 = index_q2(model, env, quat = false)
	iγ1 = index_γ1(model, env, quat = false)
	ib1 = index_b1(model, env, quat = false)
	iψ1 = index_ψ1(model, env, quat = false)
	is1 = index_s1(model, env, quat = false)
	iη1 = index_η1(model, env, quat = false)
	is2 = index_s2(model, env, quat = false)

	rz[idyn, iq2] = model.dyn.∂q2(h, q0, q1, u1, w1, λ1, q2) + s.con.dJ(λ1, q2)
	rz[idyn, [iq2; iγ1; ib1]] += s.con.dλ1(q2)' * s.con.dcf(γ1, b1, q2, k)

	# Maximum dissipation
	rz[imdp, [iψ1; iη1]] = s.con.mdψη(vT, ψ1, η1)
	rz[imdp, iq2] += s.con.mdvs(vT, ψ1, η1) * s.con.vsq2(q1, q2, k, h)

	# Other constraints
	s.con.rcz(view(rz, collect([iimp; ifri; ibimp; ibmdp; ibfri]), :), z, θ)
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

	idyn = index_dyn(model, env, quat = false)
	iimp = index_imp(model, env, quat = false)
	imdp = index_mdp(model, env, quat = false)
	ifri = index_fri(model, env, quat = false)
	ibimp = index_bimp(model, env, quat = false)
	ibmdp = index_bmdp(model, env, quat = false)
	ibfri = index_bfri(model, env, quat = false)

	iq0 = index_q0(model)
	iq1 = index_q1(model)
	iu1 = index_u1(model)
	iw1 = index_w1(model)
	iμ  = index_μ(model)
	ih  = index_h(model)

	# Dynamics
	idx = collect([iq0; iq1; iu1; iw1; ih])
	rθ[idyn, idx] = model.dyn.dθ(h, q0, q1, u1, w1, λ1, q2)

	# Maximum dissipation
	idx = collect([iq1; ih])
	rθ[imdp, idx] = s.con.vsq1h(q1, q2, k, h)

	# Other constraints
	s.con.rcθ(view(rθ, collect([iimp; ifri; ibimp; ibmdp; ibfri]), :), z, θ)
end
