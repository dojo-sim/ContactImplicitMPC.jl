@testset "set_traj!" begin
    # Test set_traj!
	T = Float64
	κ = 1e-4
	model = ContactControl.get_model("quadruped")
	ref_traj0 = ContactControl.get_trajectory("quadruped", "gait1")
	# time
	h = ref_traj0.h
	H = ref_traj0.H

	# initial conditions
	q0 = SVector{model.dim.q}(ref_traj0.q[1])
	q1 = SVector{model.dim.q}(ref_traj0.q[2])

	function ContactControl.z_initialize!(z, model::Quadruped, q1)
		nq = model.dim.q
		z .= 1.0e-1
		z[1:nq] = q1
	end

	sim0 = ContactControl.simulator(model, q0, q1, h, H,
		p = open_loop_policy(SVector{model.dim.u}.(ref_traj0.u), h),
		r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
		rz = model.spa.rz_sp,
		rθ = model.spa.rθ_sp,
		ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
		sim_opts = ContactControl.SimulatorOptions(warmstart = true))
	ContactControl.simulate!(sim0; verbose = false)

	ref_traj0 = deepcopy(sim0.traj)
	ref_traj1 = deepcopy(ref_traj0)
	for t = 1:H+2
		ref_traj1.q[t] .+= 1e-0*rand(model.dim.q)
	end

	impl0 = ContactControl.ImplicitTraj(H, model)
	impl1 = ContactControl.ImplicitTraj(H, model)
	ContactControl.linearization!(model, ref_traj0, impl0, ref_traj0.κ)
	ContactControl.linearization!(model, ref_traj0, impl1, ref_traj0.κ)

	ContactControl.implicit_dynamics!(model, ref_traj0, impl0)
	@test mean(norm.(impl0.d, 2)) < 5e-3
	ContactControl.implicit_dynamics!(model, ref_traj1, impl1)
	@test mean(norm.(impl1.d, 2)) > 5e-1



	# Check that we can optimize z with the residual function r and rz
	# Verify that we get the ~same results using r_approx and rz_approx if the linearization was done about the solution.
	nq = model.dim.q
	nz = ContactControl.num_var(model)
	α = 5e-2

	ref_traj0 = deepcopy(sim0.traj)
	ref_traj1 = deepcopy(ref_traj0)
	ref_traj1.q[3] .+= α*ones(model.dim.q)
	ref_traj1.z[1][1:nq] .= ref_traj1.q[3]

	impl0 = ContactControl.ImplicitTraj(H, model)
	ContactControl.linearization!(model, ref_traj0, impl0, ref_traj0.κ)

	z1 = ref_traj1.z[1]
	θ1 = ref_traj1.θ[1]
	@test norm(θ1 - ref_traj0.θ[1], 1) < 1e-8
	@test abs(norm(z1 - ref_traj0.z[1], 1) - nq*α) < 1e-8

	r1 = zeros(nz)
	κ = ref_traj1.κ
	rz1 = spzeros(nz,nz)
	model.res.r(r1, z1, θ1, κ[1])
	@test norm(r1) > 1.0

	function dummy_newton(z, θ, κ)
		for k = 1:400
			r = zeros(nz)
			rz = spzeros(nz,nz)
			rz = similar(model.spa.rz_sp)
			model.res.r(r, z, θ, κ)
			model.res.rz(rz, z, θ)
			Δ = - rz \ r
			z = z + 0.1*Δ
			# @show norm(r)
		end
		return z
	end

	z2 = dummy_newton(z1, θ1, κ[1])
	model.res.r(r1, z2, θ1, κ[1])
	@test norm(r1) < 1e-10

	function dummy_linear_newton(impl, z, θ, κ)
		for k = 1:400
			r = zeros(nz)
			rz = spzeros(nz,nz)
			rz = similar(model.spa.rz_sp)
			ContactControl.r_approx!(impl.lin[1], r, z, θ, κ[1])
			ContactControl.rz_approx!(impl.lin[1], rz, z, θ)
			Δ = - rz \ r
			z = z + 0.1*Δ
			# @show norm(r)
		end
		return z
	end

	z3 = dummy_linear_newton(impl0, z1, θ1, κ)
	ContactControl.r_approx!(impl0.lin[1], r1, z3, θ1, κ[1])
	@test norm(r1) < 1e-10

	# We recover the original z using r and rz
	@test norm(ref_traj0.z[1] - z2) < 1e-6
	# We recover the original z using r_approx and rz_approx
	@test norm(ref_traj0.z[1] - z3) < 1e-6
	# We recover the same solution using either methods
	@test norm(z2 - z3) < 1e-6
end

# @testset "Barrier Function" begin
# 	n = 10
# 	m = 5
# 	a = zeros(n)
# 	b = rand(SizedVector{m})
# 	v = view(a, 1:m)
# 	ContactControl.set!(v, b)
# 	@test v == b
# 	ContactControl.setminus!(v, -b)
# 	@test v == b
#
# 	a0 = rand(SizedVector{n})
# 	a1 = rand(SizedVector{n})
# 	a2 = rand(SizedVector{n})
# 	ContactControl.delta!(a0, a1, a2)
# 	@test a0 == a1 - a2
# end


@testset "Copy_traj!" begin
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

	@test traj0.q[1] == 10.0*ones(nq)
	traj0.q[1] .+= 10.0

	@test traj1.q[1] == 10.0*ones(nq)
end


@testset "Residual and Jacobian" begin
	model = ContactControl.get_model("quadruped")
	H = 2
	h = 0.1
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq+nc+nb
	nr = nq+nu+nc+nb+nd

	# Test Residual
	r0 = ContactControl.Residual(H, model.dim)
	@test length(r0.r) == H*nr

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
	j0 = ContactControl.Jacobian(H, model.dim)
	@test size(j0.j) == (H*nr,H*nr)

	# Test views
	j0.Qq2[1] .= 1.0
	j0.Qq2[2] .= 2.0
	@test all(j0.j[1:nq, 1:nq] .== 1.0)
	@test all(j0.j[nr .+ (1:nq), nr .+ (1:nq)] .== 2.0)

	# Test views
	j0.j[1:nr, 1:nr] .= -1.0
	@test all(j0.Qq2[1] .== -1.0)
	@test all(j0.Qu1[1] .== -1.0)
	@test all(j0.Qγ1[1] .== -1.0)
	@test all(j0.Qb1[1] .== -1.0)
end

@testset "jacobian!" begin
	model = ContactControl.get_model("quadruped")
	κ = 1e-4
	ref_traj0 = ContactControl.get_trajectory("quadruped", "gait1")
	ref_traj0.κ .= κ
	H = ref_traj0.H
	h = 0.1
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq+nc+nb
	nr = nq+nu+nc+nb+nd

	# Test Jacobian!
	cost0 = CostFunction(H, model.dim,
	    Qq=fill(Diagonal(1e-0*ones(SizedVector{nq})), H),
	    Qu=fill(Diagonal(1e-1*ones(SizedVector{nu})), H),
	    Qγ=fill(Diagonal(1e-2*ones(SizedVector{nc})), H),
	    Qb=fill(Diagonal(1e-3*ones(SizedVector{nb})), H),
	    )
	core0 = ContactControl.Newton(H, h, model, cost=cost0)
	impl0 = ContactControl.ImplicitTraj(H, model)
	ContactControl.linearization!(model, ref_traj0, impl0)
	ContactControl.jacobian!(model, core0, core0.j, impl0)
	spy(Matrix(core0.j.j[1:150,1:150]))

	# Test symmetry
	@test core0.j.j - core0.j.j' == spzeros(H*nr,H*nr)

	# Test cost function terms and regularization terms
	off = 0
	@test all(abs.(diag(Matrix(core0.j.j[off .+ (1:nq), off .+ (1:nq)] .- 1e-0))) .< 1e-8); off += nq
	@test all(abs.(diag(Matrix(core0.j.j[off .+ (1:nu), off .+ (1:nu)] .- 1e-1))) .< 1e-8); off += nu
	@test all(abs.(diag(Matrix(core0.j.j[off .+ (1:nc), off .+ (1:nc)] .- 1e-2))) .< 1e-8); off += nc
	@test all(abs.(diag(Matrix(core0.j.j[off .+ (1:nb), off .+ (1:nb)] .- 1e-3))) .< 1e-8); off += nb
	@test all(abs.(diag(Matrix(core0.j.j[off .+ (1:nd), off .+ (1:nd)] .+ core0.n_opts.β*impl0.lin[1].κ0[1]))) .< 1e-8); off += nd

	# Test dynamics terms
	for t = 1:H
	    @test all(diag(core0.j.IV[t]) .== -1.0)
	end
	for t = 3:H
	    @test all(core0.j.q0[t-2] .== impl0.δq0[t])
	end
	for t = 2:H
	    @test all(core0.j.q1[t-1] .== impl0.δq1[t])
	end
	for t = 1:H
	    @test all(core0.j.u1[t] .== impl0.δu1[t])
	end
end


@testset "residual!" begin
	model = ContactControl.get_model("quadruped")
	κ = 1e-4
	ref_traj0 = ContactControl.get_trajectory("quadruped", "gait1")
	ref_traj0.κ .= κ
	H = ref_traj0.H
	h = 0.1
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq+nc+nb
	nr = nq+nu+nc+nb+nd

	# Test Jacobian!
	cost0 = ContactControl.CostFunction(H, model.dim,
	    Qq=fill(Diagonal(1e-0*ones(SizedVector{nq})), H),
	    Qu=fill(Diagonal(1e-1*ones(SizedVector{nu})), H),
	    Qγ=fill(Diagonal(1e-2*ones(SizedVector{nc})), H),
	    Qb=fill(Diagonal(1e-3*ones(SizedVector{nb})), H),
	    )
	core0 = ContactControl.Newton(H, h, model, cost=cost0)
	impl0 = ContactControl.ImplicitTraj(H, model)
	ContactControl.linearization!(model, ref_traj0, impl0)
	ContactControl.implicit_dynamics!(model, ref_traj0, impl0)
	# Offset the trajectory and the dual variables to get a residual
	traj1 = deepcopy(ref_traj0)
	for t = 1:H
	    core0.ν[t] .+= 2.0
	    traj1.θ[t] .+= 0.1
	    traj1.z[t] .+= 0.1
	    traj1.u[t] .+= 0.1
	    traj1.γ[t] .+= 0.1
	    traj1.b[t] .+= 0.1
	end
	for t = 1:H+2
	    traj1.q[t] .+= 0.1
	end

	ContactControl.residual!(model, core0, core0.r, core0.ν, impl0, traj1, ref_traj0)
	@test norm(core0.Δq[1] .- 0.1, Inf) < 1e-8
	@test norm(core0.Δu[1] .- 0.1, Inf) < 1e-8
	@test norm(core0.Δγ[1] .- 0.1, Inf) < 1e-8
	@test norm(core0.Δb[1] .- 0.1, Inf) < 1e-8

	core0.r.q2[1]
	core0.cost.Qq[1]*core0.Δq[1]
	@test norm(core0.r.q2[1] .- core0.cost.Qq[1]*core0.Δq[1] - impl0.δq0[1+2]'*core0.ν[1+2] - impl0.δq1[1+1]'*core0.ν[1+1] .+ core0.ν[1][1], Inf) < 1e-8
	@test norm(core0.r.u1[1] .- core0.cost.Qu[1]*core0.Δu[1] - impl0.δu1[1]'*core0.ν[1], Inf) < 1e-8
	@test norm(core0.r.γ1[1] .- core0.cost.Qγ[1]*core0.Δγ[1] .+ core0.ν[1][1], Inf) < 1e-8
	@test norm(core0.r.b1[1] .- core0.cost.Qb[1]*core0.Δb[1] .+ core0.ν[1][1], Inf) < 1e-8
end
