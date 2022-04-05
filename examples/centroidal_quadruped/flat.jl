# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC
using LinearAlgebra

function simulate!(s::Simulator{T}; verbose=false) where T
    status = false

    N = length(s.traj.u)
    p = s.policy
    w = s.dist
    traj = s.traj

    for t = 1:N
		println("t $t $N")
        # policy
        policy_time = @elapsed traj.u[t] .= policy(p, traj, t)
		@show norm((traj.u[t] - u_ref)[[3,6,9,12]])
        s.opts.record && (s.stats.policy_time[t] = policy_time)

        # disturbances
        traj.w[t] .= disturbances(w, traj.q[t+1], t)

        # step
        status = RoboDojo.step!(s, t, verbose=verbose)
        !status && break
    end

    return status
end


# ## Simulation
s = get_simulation("centroidal_quadruped", "flat_3D_lc", "flat")
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
	# joinpath(module_dir(), "src/dynamics/centroidal_quadruped/gaits/inplace_trot_v4.jld2"),
    joinpath(module_dir(), "src/dynamics/centroidal_quadruped/gaits/stand_euler_v0.jld2"),
    load_type = :split_traj_alt));

H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 30
h_sim = h / N_sample
H_sim = 1000
κ_mpc = 2.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    v = [Diagonal(1e-3 * [[1,1,1]; 0e-1*[1,1,1]; fill([1,1,1], 4)...]) for t = 1:H_mpc],
    q = [Diagonal(3e-1 * [[3,3,10]; 1e0*[1,1,1]; 0.3fill([1,1,10], 4)...]) for t = 1:H_mpc],
    u = [Diagonal(3e-3 * vcat(fill([1,1,1], 4)...)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc]);

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
	ip_opts = InteriorPointOptions(
					undercut = 5.0,
					κ_tol = κ_mpc,
					r_tol = 1.0e-4, # TODO maybe relax this parameter
					diff_sol = true,
					solver = :empty_solver,
					max_time = 1e5),
    n_opts = NewtonOptions(
        r_tol = 3e-4,
		solver=:lu_solver,
		# solver=:ldl_solver,
        max_iter = 5),
    mpc_opts = CIMPCOptions());

# ## Disturbances
idx_d1 = 20
idx_d2 = idx_d1 + 350
idx_d3 = idx_d2 + 350
idx_d4 = idx_d3 + 350
idx_d5 = idx_d4 + 350
idx_d6 = idx_d5 + 350
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5, idx_d6]
impulses = [[0,0,10.0], [0,0,-10.0], [0,2,0.0], [0,-2,0.0], [2,0,0.0], [-2,0,0.0]]
d = impulse_disturbances(impulses, idx)

w = [[0.05,0.0,0.0] for i=1:H_sim/N_sample]
w = [[0,0.07,0.0] for i=1:H_sim/N_sample]
w = [[0,0.0,0.0] for i=1:H_sim/N_sample]
d = open_loop_disturbances(w, N_sample)

# ## Initial conditions
q1_sim, v1_sim = initial_conditions(ref_traj);

# ## Simulator

u_ref = deepcopy(ref_traj.u[1])
Δu = 0.03
u_ref[3] -= Δu
u_ref[6] += Δu
u_ref[9] -= Δu
u_ref[12] += Δu
# p = open_loop_policy([u_ref for i=1:1000], N_sample=N_sample)
sim = simulator(s, H_sim, h=h_sim, dist=d);
# sim = simulator(s, H_sim, h=h_sim, policy=p, dist=d);

using BenchmarkTools
# ## Simulate
q1_sim0 = deepcopy(q1_sim)
# Δx = 0.0
# q1_sim0[1] += Δx
# q1_sim0[7] += Δx
# q1_sim0[10] += Δx
# q1_sim0[13] += Δx
# q1_sim0[16] += Δx
# Δθ = 0.05
# q1_sim0[4] += Δθ
Δz = 0.1
q1_sim0[3] += Δz
# RoboDojo.simulate!(sim, q1_sim, v1_sim)
RoboDojo.simulate!(sim, q1_sim0, v1_sim)

# # ## Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize!(vis, model, sim.traj.q; Δt=h_sim)
# anim = visualize_meshrobot!(vis, model, sim.traj, h=h_sim * 10, sample=10);

# # ## Timing result
# # Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio
