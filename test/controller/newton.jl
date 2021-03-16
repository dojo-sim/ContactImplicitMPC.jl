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
    @test target.q[1][1] == 201.0
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

	sim0 = ContactControl.simulator2(model, q0, q1, h, H,
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
	mean(norm.(impl0.d, 2)) < 5e-3
	ContactControl.implicit_dynamics!(model, ref_traj1, impl1)
	mean(norm.(impl1.d, 2)) > 5e-1
end
