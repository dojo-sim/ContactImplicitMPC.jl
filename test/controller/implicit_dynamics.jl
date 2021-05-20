@testset "Controller: Implicit Dynamics" begin
	T = Float64
	s = get_simulation("particle", "flat_3D_lc", "flat")

	model = s.model
	env = s.env

	nq = model.dim.q
	nu = model.dim.u
	H = 10
	h = 0.0001
	κ = 1e-6
	nq = model.dim.q

	q0 = SVector{nq,T}([0.0, 0.0, 1.0])
	q1 = SVector{nq,T}([0.0, 0.0, 1.0])

	ip_opts = ContactControl.InteriorPointOptions(
		κ_init = 1.0, κ_tol = 1.0e-8, r_tol = 1.0e-10)

	sim = ContactControl.simulator(s, q0, q1, h, H;
		p = ContactControl.no_policy(model),
		d = ContactControl.no_disturbances(model),
		ip_opts = ip_opts,
		sim_opts = ContactControl.SimulatorOptions{T}())

	ContactControl.simulate!(sim, verbose = false)
	ref_traj = deepcopy(sim.traj)
	traj0 = deepcopy(ref_traj)

	im_traj = ContactControl.ImplicitTraj(ref_traj, s)
	ContactControl.implicit_dynamics!(im_traj, s, traj0)

	# Implicit dynamics contraint violation is ≈ 0
	for t = 1:H
		@test norm(im_traj.d[t], Inf) < 1.0e-5
	end

	# # We can use Newton's method to correct the trajectory and make it dynamically feasible
	# traj1 = deepcopy(ref_traj)
	# traj1.q[1] .+= 0.3
	# traj1.q[2] .+= 0.3
	# traj1.u[1] .+= 0.3
	# ContactControl.update_θ!(traj1, 1)
	# ContactControl.update_θ!(traj1, 2)
	# ContactControl.implicit_dynamics!(im_traj, s, traj1, κ = κ)
	#
	# @test norm(im_traj.d[1]) > 1.0e-5
	#
	# for k = 1:100
	# 	ContactControl.implicit_dynamics!(im_traj, s, traj1, κ = κ)
	# 	Δθ = - Matrix([im_traj.δq0[1] im_traj.δq1[1] im_traj.δu1[1]]) \ im_traj.d[1]
	#
	# 	off = 0
	# 	traj1.q[1] .+= 0.2 * Δθ[off .+ (1:nq)]; off += nq
	# 	traj1.q[2] .+= 0.2 * Δθ[off .+ (1:nq)]; off += nq
	# 	traj1.u[1] .+= 0.2 * Δθ[off .+ (1:nu)]; off += nu
	#
	# 	ContactControl.update_θ!(traj1, 1)
	# 	ContactControl.update_θ!(traj1, 2)
	# end
	#
	# ContactControl.implicit_dynamics!(im_traj, s, traj1, κ = κ)
	# @test norm(im_traj.d[1]) < 1.0e-5
end
