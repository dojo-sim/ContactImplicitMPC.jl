"""
    Feedback linearization policy for Biped5
"""

@with_kw mutable struct FeedbackLin12Options
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct FeedbackLin12 <: Policy
	kp::AbstractVector
	kd::AbstractVector
	pinnedmodel::PinnedBiped513
	opts::FeedbackLin12Options
end

function feedback_lin_policy(pinnedmodel::PinnedBiped513;
		kp = ones(4),
		kv = ones(4),
		mpc_opts = FeedbackLin12Options(),
		)

	h = ref_traj.h
	u  = zeros(model.dim.u)
	q0 = zeros(model.dim.q)
	q1 = zeros(model.dim.q)
	contact = false
	FeedbackLin12(kp, kv, pinnedmodel, mpc_opts)
end

function policy(p::FeedbackLin12, x, traj, t)
	h = traj.h
	q0 = copy(x)
	q1 = copy(traj.q[t+1])
	q  = (q0+q1)/2
	qd = (q1-q0)/h
	u_calf1, u_thigh1, u_thigh2, u_calf2 = [0.0, feedback_linearized_policy(p.pinnedmodel, [q; qd])]
	u_torso = 0.0 #TODO need fix

	u = [u_torso, u_thigh1, u_calf1, u_thigh2, u_calf2]
    return u
end
