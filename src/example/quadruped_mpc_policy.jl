# get model
model = get_model("quadruped")

# get trajectory
ref_traj = get_trajectory("quadruped", "gait2")
ref_traj_copy = deepcopy(ref_traj)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 2
H_mpc = 10
h_sim = h / N_sample
H_sim = 250

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(1e-2 * [0.02; 0.02; 1.0; 0.15*ones(model.dim.q - 3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        solver = :ldl_solver,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions())

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)

qq = []
for q in ref_traj_copy.q
    for i = 1:N_sample
        push!(qq, q)
    end
end

plot(hcat(qq...)[1:model.dim.q, 1:100]',
    label = "", color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:model.dim.q, 1:100]',
    label = "", color = :cyan, width = 1.0, legend = :topleft)

include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
vis = Visualizer()
render(vis)
anim = visualize_robot!(vis, model, sim.traj)
anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)
