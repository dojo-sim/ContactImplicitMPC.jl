@testset "Newton" begin
	T = Float64
	model = ContactControl.get_model("particle")
	H = 10
	h = 0.03
	κ = 1e-3
	nq = model.dim.q

	q0 = SVector{nq,T}([0.0, 0.0, 0.2])
	q1 = SVector{nq,T}([0.1, 0.1, 0.2])
	ip_opts = ContactControl.InteriorPointOptions(κ_init=κ, κ_tol=κ*2, r_tol=1e-5)
	sim0 = ContactControl.simulator2_base(model, q0, q1, h, H;
	    u = [@SVector zeros(model.dim.u) for t = 1:H],
	    w = [@SVector zeros(model.dim.w) for t = 1:H],
	    ip_opts = ip_opts,
	    sim_opts = ContactControl.SimulatorOptions{T}())

	ContactControl.simulate!(sim0; verbose = false)
	ref_traj0 = deepcopy(sim0.traj)
	traj0 = deepcopy(ref_traj0)

	impl0 = ContactControl.ImplicitTraj(H, model)
	ContactControl.linearization!(model, ref_traj0, impl0)
	ContactControl.implicit_dynamics!(model, traj0, impl0, κ=κ)

	# Implicit dynamics contraint violation is ≈ 0
	for k = 1:H
		@test norm(impl0.d[k], Inf) < 1e-3
	end
end
