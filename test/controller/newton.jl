@testset "Newton: update_traj!" begin
    # Test set_traj!
	T = Float64
	κ = 1.0e-4
	model = ContactControl.get_model("quadruped")
	model.μ_world = 0.5

	ref_traj = deepcopy(ContactControl.get_trajectory("quadruped", "gait2", load_type = :split_traj_alt))
	ContactControl.update_friction_coefficient!(ref_traj, model)

	# time
	h = ref_traj.h
	H = ref_traj.H

	# initial conditions
	q0 = SVector{model.dim.q}(ref_traj.q[1])
	q1 = SVector{model.dim.q}(ref_traj.q[2])

	function ContactControl.z_initialize!(z, model::ContactControl.Quadruped, q1)
		nq = model.dim.q
		z .= 1.0e-1
		z[1:nq] = q1
	end

	sim = ContactControl.simulator(model, q0, q1, h, H,
		p = ContactControl.open_loop_policy(SVector{model.dim.u}.(ref_traj.u)),
		ip_opts = ContactControl.InteriorPointOptions(
			r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
		sim_opts = ContactControl.SimulatorOptions(warmstart = true))
	ContactControl.simulate!(sim; verbose = false)

	ref_traj = deepcopy(sim.traj)
	ref_traj1 = deepcopy(ref_traj)
	for t = 1:H+2
		ref_traj1.q[t] .+= 1.0e-0 * rand(model.dim.q)
	end

	im_traj0 = ContactControl.ImplicitTraj(ref_traj, model)
	im_traj1 = ContactControl.ImplicitTraj(ref_traj1, model)

	ContactControl.implicit_dynamics!(im_traj0, model, ref_traj)
	@test maximum(norm.(im_traj0.d, 2)) < 5.0e-3
	ContactControl.implicit_dynamics!(im_traj1, model, ref_traj1)
	@test maximum(norm.(im_traj1.d, 2)) > 5.0e-1
	#
	# # Check that we can optimize z with the residual function r and rz
	# # Verify that we get the ~same results using r_approx and rz_approx if the linearization was done about the solution.
	# nq = model.dim.q
	# nz = ContactControl.num_var(model)
	# α = 5.0e-2
	#
	# ref_traj = deepcopy(sim.traj)
	# ref_traj1 = deepcopy(ref_traj)
	# ref_traj1.q[3] .+= α * ones(model.dim.q)
	# ref_traj1.z[1][1:nq] .= ref_traj1.q[3]
	#
	# im_traj0 = ContactControl.ImplicitTraj(ref_traj, model)
	#
	# z1 = ref_traj1.z[1]
	# θ1 = ref_traj1.θ[1]
	# @test LinearAlgebra.norm(θ1 - ref_traj.θ[1], 1) < 1.0e-8
	# @test abs(norm(z1 - ref_traj.z[1], 1) - nq * α) < 1.0e-8
	#
	# r1 = zeros(nz)
	# κ = ref_traj1.κ
	# rz1 = spzeros(nz,nz)
	# model.res.r!(r1, z1, θ1, κ[1])
	# @test norm(r1) > 1.0
	#
	# #TODO: add tests
	#
	# # function dummy_newton(z, θ, κ)
	# # 	for k = 1:400
	# # 		r = zeros(nz)
	# # 		rz = spzeros(nz,nz)
	# # 		rz = similar(model.spa.rz_sp)
	# # 		model.res.r!(r, z, θ, κ, nothing)
	# # 		model.res.rz!(rz, z, θ, nothing)
	# # 		Δ = - rz \ r
	# # 		z = z + 0.1 * Δ
	# # 		# @show norm(r)
	# # 	end
	# # 	return z
	# # end
	# #
	# # z2 = dummy_newton(z1, θ1, κ[1])
	# # model.res.r!(r1, z2, θ1, κ[1], nothing)
	# # @test norm(r1) < 1.0e-10
	# #
	# # function dummy_linear_newton(im_traj, z, θ, κ)
	# # 	for k = 1:400
	# # 		r = zeros(nz)
	# # 		rz = spzeros(nz, nz)
	# # 		rz = similar(model.spa.rz_sp)
	# # 		model.linearized.r!(im_traj.lin[1], r, z, θ, κ[1])
	# # 		model.linearized.rz!(im_traj.lin[1], rz, z, θ)
	# # 		Δ = - rz \ r
	# # 		z = z + 0.1 * Δ
	# # 		# @show norm(r)
	# # 	end
	# # 	return z
	# # end
	# #
	# # z3 = dummy_linear_newton(im_traj0, z1, θ1, κ)
	# # model.linearized.r!(im_traj0.lin[1], r1, z3, θ1, κ[1])
	# #
	# # @test norm(r1) < 1.0e-10
	# #
	# # # We recover the original z using r and rz
	# # @test norm(ref_traj.z[1] - z2) < 1.0e-6
	# #
	# # # We recover the original z using r_approx and rz_approx
	# # @test norm(ref_traj.z[1] - z3) < 1.0e-6
	# #
	# # # We recover the same solution using either methods
	# # @test norm(z2 - z3) < 1.0e-6
end

@testset "Newton: copy_traj!" begin
	# test copy_traj!
	model = ContactControl.get_model("quadruped")
	H = 59
	h = 0.1
	nq = model.dim.q

	traj0 = ContactControl.contact_trajectory(H, h, model)
	traj1 = ContactControl.contact_trajectory(H, h, model)
	traj1.q[1] .+= 10.0

	traj1.q
	ContactControl.copy_traj!(traj0, traj1, 1)

	@test traj0.q[1] == 10.0 * ones(nq)

	traj0.q[1] .+= 10.0
	@test traj1.q[1] == 10.0 * ones(nq)
end

@testset "Newton: Residual and Jacobian" begin
	model = ContactControl.get_model("quadruped")
	H = 2
	h = 0.1
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq + nc + nb
	nr = nq + nu + nc + nb + nd

	# Test Residual
	r0 = ContactControl.NewtonResidual(H, model.dim)
	@test length(r0.r) == H * nr

	# Test views
	r0.q2[1] .= 1.0
	r0.q2[2] .= 2.0
	@test all(r0.r[1:nq] .== 1.0)
	@test all(r0.r[nr .+ (1:nq)] .== 2.0)

	# Test views
	r0.r[1:nr] .= -1.0
	@test all(r0.q2[1] .== -1.0)
	@test all(r0.u1[1] .== -1.0)
	@test all(r0.γ1[1] .== -1.0)
	@test all(r0.b1[1] .== -1.0)

	# Test Jacobian
	j0 = ContactControl.NewtonJacobian(H, model.dim)
	@test size(j0.R) == (H * nr, H * nr)

	# Test views
	j0.obj_q2[1] .= 1.0
	j0.obj_q2[2] .= 2.0
	@test all(j0.R[1:nq, 1:nq] .== 1.0)
	@test all(j0.R[nr .+ (1:nq), nr .+ (1:nq)] .== 2.0)

	# Test views
	j0.R[1:nr, 1:nr] .= -1.0
	@test all(j0.obj_q2[1] .== -1.0)
	@test all(j0.obj_u1[1] .== -1.0)
	@test all(j0.obj_γ1[1] .== -1.0)
	@test all(j0.obj_b1[1] .== -1.0)
end

@testset "Newton: jacobian!" begin
	model = ContactControl.get_model("quadruped")
	κ = 1.0e-4

	ref_traj = get_trajectory("quadruped", "gait2", load_type = :split_traj_alt)
	ref_traj.κ .= κ
	H = ref_traj.H
	h = 0.1

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq + nc + nb
	nr = nq + nu + nc + nb + nd

	# Test Jacobian!
	obj = TrackingObjective(H, model.dim,
	    q = [Diagonal(1.0 * ones(nq)) for t = 1:H],
	    u = [Diagonal(1.0e-1 * ones(nu)) for t = 1:H],
	    γ = [Diagonal(1.0e-2 * ones(nc)) for t = 1:H],
	    b = [Diagonal(1.0e-3 * ones(nb)) for t = 1:H])
	im_traj0 = ContactControl.ImplicitTraj(ref_traj, model)
	core = ContactControl.Newton(H, h, model, ref_traj, im_traj0, obj = obj)
	ContactControl.jacobian!(core.jac, im_traj0, obj, ref_traj.H, core.β)

	# Test symmetry
	@test core.jac.R - core.jac.R' == spzeros(H * nr, H * nr)

	# Test obj function terms and regularization terms
	off = 0
	@test all(abs.(diag(Matrix(core.jac.R[off .+ (1:nq), off .+ (1:nq)] .- 1e-0))) .< 1e-8); off += nq
	@test all(abs.(diag(Matrix(core.jac.R[off .+ (1:nu), off .+ (1:nu)] .- 1e-1))) .< 1e-8); off += nu
	@test all(abs.(diag(Matrix(core.jac.R[off .+ (1:nc), off .+ (1:nc)] .- 1e-2))) .< 1e-8); off += nc
	@test all(abs.(diag(Matrix(core.jac.R[off .+ (1:nb), off .+ (1:nb)] .- 1e-3))) .< 1e-8); off += nb
	@test all(abs.(diag(Matrix(core.jac.R[off .+ (1:nd), off .+ (1:nd)] .+ core.β * im_traj0.lin[1].κ[1]))) .< 1e-8); off += nd

	# Test dynamics terms
	for t = 1:H
	    @test all(core.jac.IV[t] .== -1.0)
	end

	for t = 3:H
	    @test all(core.jac.q0[t-2] .== im_traj0.δq0[t])
	end

	for t = 2:H
	    @test all(core.jac.q1[t-1] .== im_traj0.δq1[t])
	end

	for t = 1:H
	    @test all(core.jac.u1[t] .== im_traj0.δu1[t])
	end
end

@testset "Newton: residual!" begin
	model = ContactControl.get_model("quadruped")
	κ = 1.0e-4
	ref_traj = ContactControl.get_trajectory("quadruped", "gait0")
	ref_traj.κ .= κ
	H = ref_traj.H
	h = 0.1

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq + nc + nb
	nr = nq + nu + nc + nb + nd

	# Test Jacobian!
	obj = TrackingObjective(H, model.dim,
	    q = [Diagonal(1.0 * ones(nq)) for t = 1:H],
	    u = [Diagonal(1.0 * ones(nu)) for t = 1:H],
	    γ = [Diagonal(1.0e-6 * ones(nc)) for t = 1:H],
	    b = [Diagonal(1.0e-6 * ones(nb)) for t = 1:H])
	im_traj0 = ContactControl.ImplicitTraj(ref_traj, model)
	core = ContactControl.Newton(H, h, model, ref_traj, im_traj0, obj = obj)
	ContactControl.implicit_dynamics!(im_traj0, model, ref_traj)

	# Offset the trajectory and the dual variables to get a residual
	traj1 = deepcopy(ref_traj)

	for t = 1:H
	    core.ν[t] .+= 2.0
	    traj1.θ[t] .+= 0.1
	    traj1.z[t] .+= 0.1
	    traj1.u[t] .+= 0.1
	    traj1.γ[t] .+= 0.1
	    traj1.b[t] .+= 0.1
	end

	for t = 1:H+2
	    traj1.q[t] .+= 0.1
	end

	ContactControl.residual!(core.res, model, core, core.ν, im_traj0, traj1, ref_traj)

	@test norm(core.Δq[1] .- 0.1, Inf) < 1.0e-8
	@test norm(core.Δu[1] .- 0.1, Inf) < 1.0e-8
	@test norm(core.Δγ[1] .- 0.1, Inf) < 1.0e-8
	@test norm(core.Δb[1] .- 0.1, Inf) < 1.0e-8

	@test (norm(core.res.q2[1] .- core.obj.q[1] * core.Δq[1]
		- im_traj0.δq0[1+2]' * core.ν[1+2] - im_traj0.δq1[1+1]' * core.ν[1+1]
		.+ core.ν[1][1], Inf) < 1.0e-8)
	@test (norm(core.res.u1[1] .- core.obj.u[1] * core.Δu[1]
		- im_traj0.δu1[1]' * core.ν[1], Inf) < 1.0e-8)
	@test (norm(core.res.γ1[1] .- core.obj.γ[1] * core.Δγ[1]
		.+ core.ν[1][1], Inf) < 1.0e-8)
	@test (norm(core.res.b1[1] .- core.obj.b[1] * core.Δb[1]
		.+ core.ν[1][1], Inf) < 1.0e-8)
end
