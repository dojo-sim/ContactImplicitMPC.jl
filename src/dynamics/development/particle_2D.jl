model = get_model("particle_2D")
num_var(model)
# time
h = 0.01
T = 200

# initial conditions
q0 = @SVector [0.0, 1.0]
q1 = @SVector [0.0, 1.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
status = ContactControl.simulate!(sim)

plot(hcat(sim.q[1:1:T]...)', label = ["x" "z"])
