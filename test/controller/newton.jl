@testset "Newton: copy_traj!" begin
	# test copy_traj!
	s = get_simulation("quadruped", "flat_2D_lc", "flat")
	model = s.model
	env = s.env

	H = 59
	h = 0.1
	nq = model.dim.q

	traj0 = ContactControl.contact_trajectory(model, env, H, h)
	traj1 = ContactControl.contact_trajectory(model, env, H, h)
	traj1.q[1] .+= 10.0

	ContactControl.copy_traj!(traj0, traj1, 1)

	@test traj0.q[1] == 10.0 * ones(nq)

	traj0.q[1] .+= 10.0
	@test traj1.q[1] == 10.0 * ones(nq)
end

@testset "Newton: Residual and Jacobian" begin
	s = get_simulation("quadruped", "flat_2D_lc", "flat")
	model = s.model
	env = s.env

	H = 2
	h = 0.1
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nd = nq + nc + nb
	nr = nq + nu + nc + nb + nd

	# Test Residual
	r0 = ContactControl.NewtonResidual(model, env, H)
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
	j0 = ContactControl.NewtonJacobian(model, env, H)
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
	s = get_simulation("quadruped", "flat_2D_lc", "flat")
	model = s.model
	env = s.env
	s.model.μ_world = 0.5

	ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
		joinpath(pwd(), "src/dynamics/quadruped/gaits/gait2.jld2"),
		load_type = :split_traj_alt))
	ContactControl.update_friction_coefficient!(ref_traj, s.model, s.env)

	κ = 1.0e-4
	ref_traj.κ .= κ
	H = ref_traj.H
	h = 0.1

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nd = nq + nc + nb
	nr = nq + nu + nc + nb + nd

	# Test Jacobian!
	obj = TrackingObjective(s.model, s.env, H,
	    q = [Diagonal(1.0 * ones(nq)) for t = 1:H],
	    u = [Diagonal(1.0e-1 * ones(nu)) for t = 1:H],
	    γ = [Diagonal(1.0e-2 * ones(nc)) for t = 1:H],
	    b = [Diagonal(1.0e-3 * ones(nb)) for t = 1:H])

	im_traj0 = ContactControl.ImplicitTraj(ref_traj, s)
	core = ContactControl.Newton(s, H, h, ref_traj, im_traj0, obj = obj)
	ContactControl.initialize_jacobian!(core.jac, obj, ref_traj.H)
	ContactControl.update_jacobian!(core.jac, im_traj0, obj, ref_traj.H, core.β)

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
	s = get_simulation("quadruped", "flat_2D_lc", "flat")
	model = s.model
	env = s.env
	s.model.μ_world = 0.5

	ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
		joinpath(pwd(), "src/dynamics/quadruped/gaits/gait2.jld2"),
		load_type = :split_traj_alt))
	ContactControl.update_friction_coefficient!(ref_traj, s.model, s.env)

	κ = 1.0e-4
	ref_traj.κ .= κ
	H = ref_traj.H
	h = 0.1

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nd = nq + nc + nb
	nr = nq + nu + nc + nb + nd

	# Test Jacobian!
	obj = TrackingObjective(model, env, H,
	    q = [Diagonal(1.0 * ones(nq)) for t = 1:H],
	    u = [Diagonal(1.0 * ones(nu)) for t = 1:H],
	    γ = [Diagonal(1.0e-6 * ones(nc)) for t = 1:H],
	    b = [Diagonal(1.0e-6 * ones(nb)) for t = 1:H])

	im_traj0 = ContactControl.ImplicitTraj(ref_traj, s)
	core = ContactControl.Newton(s, H, h, ref_traj, im_traj0, obj = obj)
	ContactControl.implicit_dynamics!(im_traj0, s, ref_traj)

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

	ContactControl.residual!(core.res, core, core.ν, im_traj0, traj1, ref_traj)

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

@testset "Newton: newton_solve!" begin

end
