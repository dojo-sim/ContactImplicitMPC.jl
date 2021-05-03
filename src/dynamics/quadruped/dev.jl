include(joinpath(@__DIR__, "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
render(vis)

# get hopper model
# model_sim = get_model("quadruped", surf="sinusoidal")
model = get_model("quadruped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("quadruped", "gait1", load_type=:split_traj_alt, model=model)


plot_surface!(vis, model.env, ylims=[-0.5, 0.3])
anim = visualize_meshrobot!(vis, model, ref_traj.q)
visualize_robot!(vis, model, ref_traj.q, anim=anim)

# Test visualizer
build_robot!(vis, model, name=:quad0f)
build_robot!(vis, model, name=:quad0b)
mvis = build_meshrobot!(vis, model, name=:mesh0)

t = 30
q = ref_traj.q[t] + [0,0,0.0,0,0,0,0,0,0,0,0]
q = ref_traj.q[t] + [0,0,pi/1,0,0,0,0,0,0,0,0]
set_robot!(vis, model, q, name=:quad0f, offset=0.00)
set_robot!(vis, model, q, name=:quad0b, offset=0.264)
set_meshrobot!(vis, mvis, model, q, name=:mesh0)
