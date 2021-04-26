include(joinpath(@__DIR__, "../..", "dynamics", "biped5", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
model_sim = get_model("biped5", surf="flat")
model_sim.μ_world = 1.0
# pinnedmodel = get_model("pinnedbiped5", surf="flat")
pinnedmodel = include(joinpath(@__DIR__, "model.jl"))
nq = model_sim.dim.q
nu = model_sim.dim.u
nc = model_sim.dim.c
nb = model_sim.dim.b


ref_traj = get_trajectory("biped5", "gait2", load_type=:split_traj)
# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 500

# Tstep = H*h/2
# Int(floor(Tstep/h_sim))


H_mpc = 10
im_traj = ImplicitTraj(ref_traj, model_sim)
obj = TrackingCost(H_mpc, model.dim,
    q = [Diagonal(1e-2 * [1.5, 1.5, 1.5, .15, .15]) for t = 1:H_mpc],
    u = [Diagonal(10e-1 * [1; 1; 1; ones(model_sim.dim.u-3)]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model_sim.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model_sim.dim.b)) for t = 1:H_mpc])

newton = Newton(H_mpc, h, model_sim,
    ref_traj, im_traj;
    obj = obj,
    opts = NewtonOptions())


plot(hcat([p_func(pinnedmodel, s) for s=-0.1:0.01:1.1]...)')
plot(hcat([pd_func(pinnedmodel, s) for s=-0.1:0.01:1.1]...)')
plot(hcat([pdd_func(pinnedmodel, s) for s=-0.1:0.01:1.1]...)')




# q_inip=[-0.228, 0.228, -0.1, -0.1, -0.3]
# q_midp=[-0.27,  0.2,   -0.1,  0.4, -0.4]
# q_endp=[-0.3, -0.1, -0.1, 0.228, -0.228]

# q_ini=[-0.40,  0.70, 0.0, 0.53, -0.584],
# q_mid=[-0.44,  0.54, 0.0, 0.66, -0.56 ],
# q_end=[-0.584, 0.53, 0.0, 0.70, -0.40 ],

q_inip = [-0.50, 0.60, -0.00, 0.35, -0.68]
q_midp = [-0.60, 0.50, -0.00, 0.65, -0.55]
q_endp = [-0.68, 0.35, -0.00, 0.60, -0.50]

REF = [q_inip, q_midp, q_endp]


q_ini = unpin_state(model, q_inip)
q_mid = unpin_state(model, q_midp)
q_end = unpin_state(model, q_endp)

build_robot!(vis, pinnedmodel, r=0.004)
set_robot!(vis, pinnedmodel, q_inip, r=0.004)
kinematics__4(pinnedmodel, q_inip, body=:calf_2)[2]
kinematics__1(pinnedmodel, q_inip, body=:calf_1)[1]
com_func(pinnedmodel, q_inip)[1]
set_robot!(vis, pinnedmodel, q_midp, r=0.004)
set_robot!(vis, pinnedmodel, q_endp, r=0.004)

build_robot!(vis, model, r=0.004, name=:ini)
build_robot!(vis, model, r=0.004, name=:mid)
build_robot!(vis, model, r=0.004, name=:end)
set_robot!(vis, model, q_ini, r=0.004, name=:ini)
set_robot!(vis, model, q_mid, r=0.004, name=:mid)
set_robot!(vis, model, q_end, r=0.004, name=:end)

q0_simp = q_inip
q1_simp = q0_simp + 2/H_sim * (q_midp - q0_simp)*0.2
q0_sim = SVector{nq}(unpin_state(model, q0_simp))
q1_sim = SVector{nq}(unpin_state(model, q1_simp))
set_robot!(vis, model, q0_sim, r=0.004)
set_robot!(vis, model, q1_sim, r=0.004)


# Disturbances
w = [zeros(nw) for t=1:Int(ceil(H_sim/N_sample))]
# w[75] += [5.0, -0.0]
d = open_loop_disturbances(w)

# Select the initial speed
p = feedback_lin_policy(pinnedmodel, kp=500*ones(4), kd=-200*ones(4))

sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    d = d,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)


plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

plot_lines!(vis, model_sim, sim.traj.q[1:10:end])
plot_surface!(vis, model_sim.env, n=200)
anim = visualize_robot!(vis, model_sim, sim.traj.q)
# anim = visualize_robot!(vis, model_sim, sim.traj, sample=5)
# anim = visualize_force!(vis, model_sim, sim.traj, anim=anim, h=h_sim, sample=5)

plot(hcat([q[[1,2,3,4,5,6,7,]] for q in Vector.(sim.traj.q)[1:end]]...)')
plot(hcat([q[[3]] for q in Vector.(sim.traj.q)[1:end]]...)')
plot(hcat([γ[[1,2]] for γ in Vector.(sim.traj.γ)[1:end]]...)')
plot(hcat([u[[1,2,3,4,5]] for u in Vector.(sim.traj.u)[1:end]]...)')


# filename = "biped_fu_policy"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
