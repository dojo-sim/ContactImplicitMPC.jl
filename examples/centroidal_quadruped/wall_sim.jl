# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC
using LinearAlgebra
using Quaternions

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

@show Threads.nthreads()

# ## Simulation
s = get_simulation("centroidal_quadruped_wall", "flat_3D_lc", "flat")
model = s.model
env = s.env

q1_sim = zeros(model.nq)
v1_sim = zeros(model.nq) 

v1_sim[1] = 10.0
v1_sim[1] = 10.0
v1_sim[7] = 10.0
v1_sim[10] = 10.0
v1_sim[13] = 10.0
v1_sim[16] = 10.0

q1_sim[3] = 1.0
q1_sim[9] = 1.0
q1_sim[12] = 1.0
q1_sim[15] = 1.0
q1_sim[18] = 1.0

# ## Simulator
sim = simulator(s, 100, h=0.1)

# ## Simulate
q1_sim0 = deepcopy(q1_sim)
RoboDojo.simulate!(sim, q1_sim0, v1_sim)

# ## Visualize
set_light!(vis)
set_floor!(vis, grid=true)
set_background!(vis)
anim = visualize!(vis, model, sim.traj.q; Δt=0.1)

sim.ip.z
sim.ip.θ



sim.traj.q[1][16]
ϕ_func(model, env, sim.traj.q[end])
