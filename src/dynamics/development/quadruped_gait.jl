# Reference trajectory
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait_1.jld2") z̄ x̄ ū q̄ τ̄ λ̄ b̄ h̄


# time
h = mean(h̄)
T = length(τ̄)

[norm(dynamics(model, h̄[t], q̄[t], q̄[t+1], h̄[t] * τ̄[t], zeros(model.dim.w), h̄[t] * λ̄[t], h̄[t] * b̄[t], q̄[t+2]), Inf) for t = 1:T]


# initial conditions
q0 = SVector{model.dim.q}(q̄[1])
q1 = SVector{model.dim.q}(q̄[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q

    z .= 3.0e-2
    z[1:nq] = q1
end

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    u = [SVector{model.dim.u}(h * u) for u in τ̄],
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = rz_sp,
    rθ = rθ_sp,
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-6, κ_tol = 1.0e-8, κ_init = 1.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

step!(sim, 1)
# simulate
status = ContactControl.simulate!(sim, verbose = true)

include(joinpath(pwd(), "src/dynamics/quadruped/visuals.jl"))
vis = Visualizer()
render(vis)
# visualize!(vis, model, q̄, Δt = h)
visualize!(vis, model, sim.q, Δt = h)
