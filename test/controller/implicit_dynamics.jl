@testset "Controller: Implicit Dynamics" begin
	T = Float64
	model = ContactControl.get_model("particle")
	nq = model.dim.q
	nu = model.dim.u
	H = 10
	h = 0.03
	κ = 1e-6
	nq = model.dim.q

	q0 = SVector{nq,T}([0.0, 0.0, 0.2])
	q1 = SVector{nq,T}([0.1, 0.1, 0.2])
	ip_opts = ContactControl.InteriorPointOptions(κ_init=κ, κ_tol=2κ, r_tol=1e-11)
	sim0 = ContactControl.simulator(model, q0, q1, h, H;
		p = ContactControl.no_policy(model),
		d = ContactControl.no_disturbances(model),
		ip_opts = ip_opts,
		sim_opts = ContactControl.SimulatorOptions{T}()
		)

	ContactControl.simulate!(sim0; verbose = false)
	ref_traj0 = deepcopy(sim0.traj)
	traj0 = deepcopy(ref_traj0)

	impl0 = ContactControl.ImplicitTraj(H, model)
	ContactControl.linearization!(model, ref_traj0, impl0)
	ContactControl.implicit_dynamics!(model, traj0, impl0, κ = κ)

	# Implicit dynamics contraint violation is ≈ 0
	for k = 1:H
		@test norm(impl0.d[k], Inf) < 1e-6
	end

	# We can use Newton's method to correct the trajectory and make it dynamically feasible
	traj1 = deepcopy(ref_traj0)
	traj1.q[1] .+= 0.3
	traj1.q[2] .+= 0.3
	traj1.u[1] .+= 0.3
	ContactControl.update_θ!(traj1, 1)
	ContactControl.update_θ!(traj1, 2)

	ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
	@test norm(impl0.d[1]) > 1e-6
	for k = 1:100
		ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
		Δθ = - Matrix([impl0.δq0[1] impl0.δq1[1] impl0.δu1[1]]) \ impl0.d[1]
		off = 0
		traj1.q[1] .+= 0.2*Δθ[off .+ (1:nq)]; off += nq
		traj1.q[2] .+= 0.2*Δθ[off .+ (1:nq)]; off += nq
		traj1.u[1] .+= 0.2*Δθ[off .+ (1:nu)]; off += nu
		ContactControl.update_θ!(traj1, 1)
		ContactControl.update_θ!(traj1, 2)
	end
	ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
	@test norm(impl0.d[1]) < 1e-6
end

#
# T = Float64
# model = ContactControl.get_model("particle")
# nq = model.dim.q
# nu = model.dim.u
# H = 10
# h = 0.03
# κ = 1e-6
# nq = model.dim.q
#
# q0 = SVector{nq,T}([0.0, 0.0, 0.2])
# q1 = SVector{nq,T}([0.1, 0.1, 0.2])
# ip_opts = ContactControl.InteriorPointOptions(κ_init=κ, κ_tol=2κ, r_tol=1e-11)
# sim0 = ContactControl.simulator(model, q0, q1, h, H;
# 	p = ContactControl.no_policy(model),
# 	d = ContactControl.no_disturbances(model),
# 	ip_opts = ip_opts,
# 	sim_opts = ContactControl.SimulatorOptions{T}()
# 	)
#
# ContactControl.simulate!(sim0; verbose = false)
# ref_traj0 = deepcopy(sim0.traj)
# traj0 = deepcopy(ref_traj0)
#
# impl0 = ContactControl.ImplicitTraj(H, model)
# ContactControl.linearization!(model, ref_traj0, impl0)
# ContactControl.implicit_dynamics!(model, traj0, impl0, κ = κ)
#
# # Implicit dynamics contraint violation is ≈ 0
# for k = 1:H
# 	@test norm(impl0.d[k], Inf) < 1e-6
# end
#
# # We can use Newton's method to correct the trajectory and make it dynamically feasible
# ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
# traj1 = deepcopy(ref_traj0)
# traj1.q[1] .+= 0.3
# traj1.q[2] .+= 0.3
# traj1.u[1] .+= 0.3
# ContactControl.update_θ!(traj1, 1)
# ContactControl.update_θ!(traj1, 2)
#
# ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
# @test norm(impl0.d[1]) > 1e-6
# for k = 1:100
# 	ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
# 	Δθ = - Matrix([impl0.δq0[1] impl0.δq1[1] impl0.δu1[1]]) \ impl0.d[1]
# 	off = 0
# 	traj1.q[1] .+= 0.2*Δθ[off .+ (1:nq)]; off += nq
# 	traj1.q[2] .+= 0.2*Δθ[off .+ (1:nq)]; off += nq
# 	traj1.u[1] .+= 0.2*Δθ[off .+ (1:nu)]; off += nu
# 	ContactControl.update_θ!(traj1, 1)
# 	ContactControl.update_θ!(traj1, 2)
# end
# ContactControl.implicit_dynamics!(model, traj1, impl0, κ = κ)
# @test norm(impl0.d[1]) < 1e-6
#
#
#
