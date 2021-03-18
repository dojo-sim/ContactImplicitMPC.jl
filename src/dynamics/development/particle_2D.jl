model = get_model("particle_2D")

# time
h = 0.01
T = 200

# initial conditions
q0 = @SVector [0.0, 1.0]
q1 = @SVector [0.0, 1.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot(hcat(sim.traj.q[1:1:T]...)', label = ["x" "z"])
