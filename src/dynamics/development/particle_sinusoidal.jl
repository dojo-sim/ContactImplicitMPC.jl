model = get_model("particle", surf = "sinusoidal")
ref_traj = ContactControl.get_trajectory("particle", "sinusoidal2",
    model_name = "particle_sinusoidal")

ref_traj.κ[1] = 1.0e-4
rep_traj = repeat_ref_traj(ref_traj, model, 2, idx_shift = (1:1))
rep_traj_copy = repeat_ref_traj(ref_traj, model, 2, idx_shift = (1:1))

T = rep_traj.H
h = rep_traj.h

# initial conditions
q0 = SVector{model.dim.q}(rep_traj.q[1])
q1 = SVector{model.dim.q}(rep_traj.q[2])

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in rep_traj.u], h, N_sample = 1),
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot(hcat(rep_traj.q...)[1:3, :]',
    label = ["x" "y" "z"], color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:3, :]',
    label = ["x" "y" "z"], color = :red, width = 1.0, legend = :topleft)

# linearized motion planning
cost = ContactControl.CostFunction(T, model.dim,
    q = [Diagonal(1.0 * ones(model.dim.q))    for t = 1:T],
    u = [Diagonal(1.0e-1 * ones(model.dim.u)) for t = 1:T],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:T],
    b = [Diagonal(1.0e-3 * ones(model.dim.b)) for t = 1:T])
solver = ContactControl.linear_motion_planning_solver(model, rep_traj, cost)
status = ContactControl.lmp!(solver)
@assert status
plot!(hcat(solver.traj.q...)[1:3, :]',
    label = ["x" "y" "z"], color = :cyan, width = 1.0, legend = :topleft)
# solver.traj.q[1] .+= [0.1; 0.1; 0.0]
# status = ContactControl.lmp!(solver)
# plot!(hcat(solver.traj.q...)[1:2, :]',
#     label = ["x" "y"], color = :magenta, width = 1.0, legend = :topleft)

tf = h * T
t_ref = range(0.0, stop = tf, length = T + 1)
length(t_ref)
H_mpc = 10
N_sample = 5
h_sim = h / N_sample
t_sim = range(0, stop = tf, length = (H * N_sample + 1))
T_sim = 300
q1_sim = SVector{model.dim.q}(copy(rep_traj.q[2]))
q0_sim = SVector{model.dim.q}(copy(q1_sim - (rep_traj.q[2] - rep_traj.q[1]) / N_sample))
(q1_sim - q0_sim) / h_sim
(rep_traj.q[2] - rep_traj.q[1]) / rep_traj.h
# sub-trajectory indices
t_idx = [collect((t - 1) .+ (1:H_mpc)) for t = 1:(H - H_mpc + 1)]

# sub-trajectories
sub_traj = [sub_ref_traj(rep_traj, model, idx) for idx in t_idx]

# create lmp solver for each MPC time step
lmp_solvers = [linear_motion_planning_solver(model, st, cost) for st in sub_traj]

# pre-solve each linearized motion planning problem
for lmp in lmp_solvers
    @show lmp!(lmp)
end

"""
    linearized model-predictive control policy
"""
struct LinearizedMPC <: Policy
   lmp
   t
   s
   q0
end

function linearized_mpc_policy(lmp, t, q0)
    LinearizedMPC(lmp, t, [i == 1 ? true : false for i = 1:length(lmp)], q0)
end

function policy(p::LinearizedMPC, x, traj, t)
    k = searchsortedlast(p.t, t)
    @show k
    if !(p.s[k] == true)
        println("re-optimizing! t = $t")
        p.lmp[k].traj.q1 .= copy(x)
        p.lmp[k].traj.q0 .= copy(p.q0)
        p.s[k] = lmp!(p.lmp[k])
        p.q0 .= copy(x)
    end
    return p.lmp[k].traj.u[1] ./ N_sample
end

p = linearized_mpc_policy(lmp_solvers, t_ref, Array(q1_sim))
#
# policy(p, nothing, sim.traj, 0.0)
# searchsortedlast(p.t, t_sim[11])

# simulator
sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, T_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

@time status = ContactControl.simulate!(sim)

qq = []
for q in rep_traj_copy.q
    for i = 1:N_sample
        push!(qq, q)
    end
end
plot(hcat(qq...)[1:3, 1:T_sim]',
    label = ["x" "y" "z"], color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:3, 1:T_sim]',
    label = ["x" "y" "z"], color = :cyan, width = 1.0, legend = :topleft)
