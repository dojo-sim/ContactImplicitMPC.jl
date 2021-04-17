"""
    Feedback linearization policy for Biped5
"""

@with_kw mutable struct FeedbackLin14Options
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct FeedbackLin14 <: Policy
	kp::AbstractVector
	kd::AbstractVector
	body::Symbol
	impact_threshold::Real
	pinnedmodel::PinnedBiped513
	opts::FeedbackLin14Options
end

function feedback_lin_policy(pinnedmodel::PinnedBiped513;
		kp = ones(4),
		kd = ones(4),
		mpc_opts = FeedbackLin14Options(),
		)

	h = ref_traj.h
	body = :toe_1
	impact_threshold = 1.0
	u  = zeros(model.dim.u)
	q0 = zeros(model.dim.q)
	q1 = zeros(model.dim.q)
	contact = false
	FeedbackLin14(kp, kd, body, impact_threshold, pinnedmodel, mpc_opts)
end

function feedback_linearized_policy(model::PinnedBiped513, x, kp, kd)
	nq = model.dim.q
	nh = 4
	q  = x[1:nq]
	qd = x[nq .+ (1:nq)]

	∇fu = ∇fu_func(model, x)
	Afu = Afu_func(model, x, ∇fu=∇fu)
	bfu = bfu_func(model, x, ∇fu=∇fu)
	h = h_func(model, q)
	hd = hd_func(model, q, qd)

	s = s_func(model, q)
	sd = sd_func(model, q, qd)
	# @show s
	# @show h

	p = p_func(model, s)
	pd = pd_func(model, s)
	pdd = pdd_func(model, s)
	# @show p

	hdd = zeros(nh)
	for i = 1:nh-1
		hdd[i] = kp[i] * (p[i] - h[i]) + kd[i]*(pd[i]*sd - hd[i]) + pdd[i]*sd^2
	end
	# #TODO need to implement the torso control
	# i = nh
	# hdd[i] = kp[i] * (p[i](s) - h[i]) + kd[i]*(pd[i](s)*sd - hd[i]) + pdd*sd^2
	# @show cond(Afu)
	@show norm(bfu)
	@show norm(hdd)
	up = Afu\(-bfu + hdd) #TODO
	return up
end


function policy(p::FeedbackLin14, x, traj, t)
	h = traj.h
	γ = traj.γ[min(traj.H, max(t-1, 1))]
	p.body = switch_foot(γ, p.body, p.impact_threshold)
	println(t, p.body)
	q0p = copy(pin_state(x, body=p.body))
	q1p = copy(pin_state(traj.q[t+1], body=p.body))
	qp  = (q0p + q1p)/2
	qdp = (q1p - q0p)/h

	# u_torso = 0.0 #TODO need fix
	up = feedback_linearized_policy(p.pinnedmodel, [qp; qdp], p.kp, p.kd)
	u = h*unpin_control(up, body=p.body)
	u = clamp.(u, -0.04, 0.04)
    return u
end


function switch_foot(γ::AbstractVector, body::Symbol, impact_threshold::Real)
	new_body = body
	if body == :toe_1
		if γ[2] > impact_threshold
			new_body = :toe_2
		end
	elseif body == :toe_2
		if γ[1] > impact_threshold
			new_body = :toe_1
		end
	end
	return new_body
end
