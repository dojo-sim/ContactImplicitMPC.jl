# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model 
env = s.env

using InteractiveUtils

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))

H = ref_traj.H
h = ref_traj.h

# ## MPC setup 
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 1000
κ_mpc = 2.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [0.1,3,1,3])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-3, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.nc)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.nc * friction_dim(env))) for t = 1:H_mpc])

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
	ip_opts = InteriorPointOptions(
					undercut = 5.0,
					κ_tol = κ_mpc,
					r_tol = 1.0e-8,
					diff_sol = true,
					solver = :empty_solver,
					max_time = 1e5),
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = CIMPCOptions(live_plotting=false))

# ## Initial conditions
q1_sim = ContactImplicitMPC.SVector{model.nq}(copy(ref_traj.q[2]))
q0_sim = ContactImplicitMPC.SVector{model.nq}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample));
v1_sim = (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / ref_traj.h

# ## Simulator
sim = Simulator(s, 10, h=h_sim)#, policy=p)

## BENCHMARK 
using InteractiveUtils


t = 1
update!(p.im_traj, p.traj, p.s, p.altitude, κ = p.κ[1]) 
@benchmark update!($p.im_traj, $p.traj, $p.s, $p.altitude, κ = $p.κ[1]) 
@code_warntype update!(p.im_traj, p.traj, p.s, p.altitude, κ = p.κ[1]) 

A = @SMatrix ones(10, 10) 

Aa = rand(10, 10)

@benchmark $A = $Aa
B = @SMatrix ones(10, 10) 
C = @SMatrix ones(10, 10) 
D = @SMatrix ones(10, 10)

M = [A B; C D]
M1 = view(M, collect(1:5), collect(1:5))
typeof(M1)
# cache = deepcopy(p.traj)
# @benchmark rotate!($p.traj, $p.traj_cache)
# @benchmark mpc_stride!($p.traj, $p.stride)

# @code_warntype rot_n_stride!(p.traj, p.traj_cache, p.stride)
# @benchmark rot_n_stride!($p.traj, $p.traj_cache, $p.stride)

# @code_warntype live_plotting(p.s.model, p.traj, sim.traj, p.newton, p.q0, sim.traj.q[t+1], t) 

policy(p, sim.traj, t)
@benchmark policy($p, $sim.traj, $t)
@code_warntype policy(p, sim.traj, t)

# @benchmark $p.ϕ .= ϕ_func($s.model, $s.env, $sim.traj.q[1])
# @benchmark s.ϕ($p.ϕ, $sim.traj.q[1])
# @benchmark ContactImplicitMPC.update_altitude!($p.altitude, $p.ϕ, $p.s,
#                  $sim.traj, $t, $p.s.model.nc, $p.N_sample,
#                  threshold = $p.opts.altitude_impact_threshold,
#                  verbose = $p.opts.altitude_verbose)

# @code_warntype ContactImplicitMPC.update_altitude!(p.altitude, p.ϕ, p.s,
#     sim.traj, t, p.s.model.nc, p.N_sample,
#     threshold = p.opts.altitude_impact_threshold,
#     verbose = p.opts.altitude_verbose)

# ## Simulate
simulate!(sim, Array(q1_sim), Array(v1_sim))

q1a = Array(q1_sim) 
v1a = Array(v1_sim)

using BenchmarkTools
@benchmark simulate!($sim, $q1a, $v1a)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_robot!(vis, model, sim.traj, h=h_sim * 5, sample=5);

# ## Timing result
# Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio

