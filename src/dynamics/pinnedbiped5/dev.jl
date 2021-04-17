include("../biped5/model.jl")
include("../biped5/visuals.jl")
include("../pinnedbiped5/model.jl")
include("../pinnedbiped5/visuals.jl")
include("../pinnedbiped5/inverted_pendulum.jl")

vis = Visualizer()
open(vis)

model = get_model("biped5")
pinnedmodel = pinnedbiped5
nq = model.dim.q
nqp = pinnedmodel.dim.q
nh = 4
build_robot!(vis, model, r=0.004)
build_robot!(vis, pinnedmodel, r=0.004)


q0 = [0.0, 0.0,-0.1,0.2,-0.3,-0.4,-0.5]
q̇0 = ones(nq)
q0p = pin_state(q0)
q̇0p = pin_state(q̇0)
set_robot!(vis, model, q0, r=0.004)
set_robot!(vis, pinnedmodel, q0p, r=0.004)


h0 = h_func(model, q0)
q0p = pin_state(q0)
h0p = h_func(pinnedmodel, q0p)
@test norm(h0 - h0p, Inf) < 1e-10

u0 = ones(nh)
x0p = [q0p; q̇0p]
f0 = f_func(pinnedmodel, x0p)
g0 = g_func(pinnedmodel, x0p)
ẋ0 = ẋ_func(pinnedmodel, x0p, u0)
∇fu0 = ∇fu_func(pinnedmodel, x0p)
@time ∇fu0 = ∇fu_func(pinnedmodel, x0p)
@time Afu0 = Afu_func(pinnedmodel, x0p)
@time bfu0 = bfu_func(pinnedmodel, x0p)

α = 0.228
qinip = [-α, α, -0.1, -0.1, -0.3]
kinematics__4(pinnedmodel, qinip, body=:calf_2, mode=:ee)[2]
set_robot!(vis, pinnedmodel, qinip, r=0.004)

qmidp = [-0.35, 0.2, -0.1, 0.3, -0.4]
kinematics__4(pinnedmodel, qmidp, body=:calf_2, mode=:ee)[2]
set_robot!(vis, pinnedmodel, qmidp, r=0.004)

qendp = [-0.3, -0.1, -0.1, α, -α]
kinematics__4(pinnedmodel, qendp, body=:calf_2, mode=:ee)[2]
set_robot!(vis, pinnedmodel, qendp, r=0.004)

θcom_ini = θcom_func(pinnedmodel, qinip)
θcom_mid = θcom_func(pinnedmodel, qmidp)
θcom_end = θcom_func(pinnedmodel, qendp)

xinip = [qinip; zeros(nqp)]
up = feedback_linearized_policy(pinnedmodel, xinip)





# ref_traj.H
# ref_traj.q[1] - ref_traj.q[end-1]
#
# function inverted_pendulum_traj(model::PinnedBiped513, traj::ContactTraj)
#     H = traj.H
#     q = ref_traj.q[1:H+1]
#     θcom_traj = zeros(H+1)
#     s_traj = zeros(H+1)
#
#     for t = 1:H+1
#         θcom_traj[t] = θcom_func(model, pin_state(q[t]))
#     end
#     θcom_ini = θcom_traj[1]
#     θcom_end = θcom_traj[end]
#     for t = 1:H+1
#         s_traj[t] = (θcom_traj[t] - θcom_ini)/(θcom_end - θcom_ini)
#     end
#     return θcom_traj, s_traj
# end
