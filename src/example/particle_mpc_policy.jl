# get model
model = get_model("particle", surf = "sinusoidal")

# get ref. trajectory
ref_traj = ContactControl.get_trajectory("particle", "sinusoidal2",
    model_name = "particle_sinusoidal")

ref_traj_copy = deepcopy(ref_traj)

# time
H = ref_traj.H
h = ref_traj.h

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[2])

# simulator
sim = ContactControl.simulator(model, q0, q1, 1.0 * h, H,
    p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u], h, N_sample = 1),
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_init = 1.0e-5, κ_tol = 1.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot(hcat(ref_traj.q...)[1:3, :]',
    label = ["x" "y" "z"], color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:3, :]',
    label = ["x" "y" "z"], color = :red, width = 1.0, legend = :topleft)

# linearized motion planning
cost = ContactControl.CostFunction(H, model.dim,
    q = [Diagonal(1.0 * ones(model.dim.q))    for t = 1:H],
    u = [Diagonal(1.0e-1 * ones(model.dim.u)) for t = 1:H],
    γ = [Diagonal(1.0e-6 * ones(model.dim.c)) for t = 1:H],
    b = [Diagonal(1.0e-6 * ones(model.dim.b)) for t = 1:H])

# Simulate linearized MPC policy
κ_mpc = 1.0e-4
tf = h * H
H_mpc = 20
N_sample = 5
h_sim = h / N_sample
H_sim = 500
q1_ref = copy(ref_traj.q[2])
q0_ref = copy(copy(ref_traj.q[1]))
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc)

# # simulator
sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-5,
        κ_tol = 1.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)

qq = []
for q in ref_traj_copy.q
    for i = 1:N_sample
        push!(qq, q)
    end
end

L = min(H_sim, length(qq))
plot(hcat(qq...)[1:3, 1:L]',
    label = ["x" "y" "z"], color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:3, 1:L]',
    label = ["x" "y" "z"], color = :cyan, width = 1.0, legend = :topleft)
