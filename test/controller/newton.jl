@testset "set_traj!" begin
    # Test set_traj!
    T = Float64
    H = 10
    h = 0.1
    model = ContactControl.get_model("quadruped")
    target = ContactControl.contact_trajectory(H, h, model)
    source = ContactControl.contact_trajectory(H, h, model)
	nq = model.dim.q # configuration
	nu = model.dim.u # control
	nc = model.dim.c # contact
	nb = model.dim.b # linear friction
	nd = nq + nc + nb # implicit dynamics constraint
	nr = nq+nu+nc+nb+nd # size of a one-time-step block
    νtarget = [-30*ones(SizedVector{nd,T}) for t=1:H]
    νsource = [+30*ones(SizedVector{nd,T}) for t=1:H]
    for t = 1:H
        source.q[t] .= +1.0
        source.u[t] .= +2.0
        source.w[t] .= +3.0
        source.γ[t] .= +4.0
        source.b[t] .= +5.0
        source.z[t] .= +6.0
        source.θ[t] .= +7.0

        target.q[t] .= -1.0
        target.u[t] .= -2.0
        target.w[t] .= -3.0
        target.γ[t] .= -4.0
        target.b[t] .= -5.0
        target.z[t] .= -6.0
        target.θ[t] .= -7.0
    end
    Δ0 = ContactControl.Residual11(H, model.dim)
    Δ0.r .+= 100*ones(nr*H)
    ContactControl.set_traj!(target, source, νtarget, νsource, Δ0, 2.0)
    source.q[1] .+= 1000.0
    source.u[1] .+= 1000.0
	@test target.q[1][1] == -1.0
	@test target.q[2][1] == -1.0
    @test target.q[3][1] == 201.0
    @test target.u[1][1] == 202.0
    @test target.γ[1][1] == 204.0
    @test target.b[1][1] == 205.0
    @test νtarget[1][1]  == 230.0
end

@testset "linearization!, implicit_dynamics!" begin
	T = Float64
	κ = 1e-4
	model = ContactControl.get_model("quadruped")
	@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b
	# time
	h = h̄
	H = length(u)

	# initial conditions
	q0 = SVector{model.dim.q}(q[1])
	q1 = SVector{model.dim.q}(q[2])

	function ContactControl.z_initialize!(z, model::Quadruped, q1)
		nq = model.dim.q
	    z .= 1.0e-1
	    z[1:nq] = q1
	end

	sim0 = ContactControl.simulator(model, q0, q1, h, H,
	    u = [SVector{model.dim.u}(h * ut) for ut in u],
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
	model.res.r(r1, z1, θ1, κ)
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

	z2 = dummy_newton(z1, θ1, κ)
	model.res.r(r1, z2, θ1, κ)
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

@testset "Barrier Function" begin
	n = 10
	m = 5
	a = zeros(n)
	b = rand(SizedVector{m})
	v = view(a, 1:m)
	ContactControl.set!(v, b)
	@test v == b
	ContactControl.setminus!(v, -b)
	@test v == b

	a0 = rand(SizedVector{n})
	a1 = rand(SizedVector{n})
	a2 = rand(SizedVector{n})
	ContactControl.delta!(a0, a1, a2)
	@test a0 == a1 - a2
end

# const ContactControl = Main


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
