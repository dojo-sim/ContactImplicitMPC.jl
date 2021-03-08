include(joinpath(pwd(), "src/dynamics/quadruped/model.jl"))

@load joinpath(res_path, "quadruped/sparse_jacobians.jld2") rz_sp rθ_sp

# time
h = 0.01
T = 100

## DROP
# initial conditions
q0 = @SVector [0.0, 0.0, 1.0]
q1 = @SVector [0.0, 0.0, 1.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = rz_sp,
    rθ = rθ_sp,
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
status = ContactControl.simulate!(sim)
@test status
@test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-6))
