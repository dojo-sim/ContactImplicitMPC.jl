# Reference trajectory
model = get_model("hopper_2D", surf = "flat")
q, u, γ, b, h = get_gait("hopper_2D", "vertical")

# time
T = length(u)

maximum([norm(ContactControl.dynamics(model,
	h, q[t], q[t+1], u[t],
	zeros(model.dim.w), γ[t], b[t], q[t+2]), Inf) for t = 1:T])

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    p = open_loop_policy([SVector{model.dim.u}(ut) for ut in u], h),
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-5, κ_init = 1.0e-4, solver = :mgs_solver),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
@time status = ContactControl.simulate!(sim, verbose = false)

include(joinpath(pwd(), "src/dynamics/hopper_2D/visuals.jl"))
vis = Visualizer()
# open(vis)
render(vis)
visualize!(vis, model, q, Δt = h)
# visualize!(vis, model, sim.traj.q, Δt = h)
