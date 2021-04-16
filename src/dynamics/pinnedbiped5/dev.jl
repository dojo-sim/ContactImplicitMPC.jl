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

function feedback_linearized_policy(model::PinnedBiped513, x)
	nq = model.dim.q
	nh = 4
	q  = x[1:nq]
	qd = x[nq .+ (1:nq)]

	∇fu = ∇fu_func(model, x)
	Afu = Afu_func(model, x, ∇fu=∇fu)
	bfu = bfu_func(model, x, ∇fu=∇fu)
	h = h_func(model, q)
	hd = hd_func(model, q, qd)

	kp = ones(nh)
	kd = ones(nh)
	s = s_func(model, q)
	sd = sd_func(model, q, qd)
	@show s
	@show h

	p = p_func(model, s)
	pd = pd_func(model, s)
	pdd = pdd_func(model, s)
	@show p

	hdd = zeros(nh)
	for i = 1:nh-1
		hdd[i] = kp[i] * (p[i] - h[i]) + kd[i]*(pd[i]*sd - hd[i]) + pdd[i]*sd^2
	end
	# #TODO need to implement the torso control
	# i = nh
	# hdd[i] = kp[i] * (p[i](s) - h[i]) + kd[i]*(pd[i](s)*sd - hd[i]) + pdd*sd^2

	up = Afu\(-bfu + hdd)
	return up
end

xinip = [qinip; zeros(nqp)]
up = feedback_linearized_policy(pinnedmodel, xinip)
