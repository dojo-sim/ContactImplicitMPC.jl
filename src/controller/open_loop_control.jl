include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper model
model = get_model("hopper_2D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

ref_traj = get_trajectory("hopper_2D", "gait_in_place", load_type=:joint_traj)
H = ref_traj.H
h = ref_traj.h
κ = 1.0e-8
