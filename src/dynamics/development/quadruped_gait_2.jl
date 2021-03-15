# Reference trajectory
model = get_model("quadruped", surf = "flat")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b

# time
h = h̄
T = length(u)

maximum([norm(dynamics(model,
	h, q[t], q[t+1], h * u[t],
	zeros(model.dim.w), h * γ[t], h * b[t], q[t+2]), Inf) for t = 1:T])

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q

    z .= 1.0e-1
    z[1:nq] = q1
end

# simulator
sim = simulator(model, q0, q1, h, T,
    u = [SVector{model.dim.u}(h * ut) for ut in u],
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-5, κ_init = 1.0e-4),
    sim_opts = SimulatorOptions(warmstart = true))

# simulate
@time status = simulate!(sim, verbose = true)

include(joinpath(pwd(), "src/dynamics/quadruped/visuals.jl"))
# vis = Visualizer()
# open(vis)
visualize!(vis, model, q, Δt = h)
visualize!(vis, model, sim.q, Δt = h)

plot(norm.(q .- sim.q))


nb
nc
transpose(sum(reshape(b[1], (Int(nb/nc), nc)), dims = 1))
