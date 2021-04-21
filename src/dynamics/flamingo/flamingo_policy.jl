"""
    Flamingo18 hopper policy for Hopper2D
"""

@with_kw mutable struct Flamingo18Options
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct Flamingo18 <: Policy
	model::Flamingo
	h_sim::Real
	u::AbstractVector
	contact::Vector{Bool}
	q0::AbstractVector
	q1::AbstractVector
	qref::AbstractVector
	xdref::Real
	opts::Flamingo18Options
end

function flamingo_policy(model::Flamingo, h_sim;
		qref=[0.0, 0.87, 0.00, 0.0, -0.0, 0.0, -0.0, -π/2, -π/2],
		xdref=0.4,
		mpc_opts = Flamingo18Options(),
		)

	u  = zeros(model.dim.u)
	q0 = zeros(model.dim.q)
	q1 = zeros(model.dim.q)
	contact = [false for i=1:model.dim.c]
	Flamingo18(model, h_sim, u, contact, copy(q0), copy(q1), qref, xdref, mpc_opts)
end

function policy(p::Flamingo19, x, traj, t)
	nc = p.model.dim.c
	il1 = [2,3,6]
	il2 = [4,5,7]

	m_flamingo = model.m_torso + model.m_thigh1 +
		model.m_thigh2 + model.m_calf1 + model.m_calf2
	# Initialization
	h = traj.h
	p.q0 = copy(traj.q[t])
	p.q1 = copy(x)

	# Detect contact
	p.contact = detect_contact(model, traj.γ[max(1,t-1)])
	@show p.contact

	# Velocities
	qd = (p.q1 - p.q0)/traj.h

	iθ = Vector(3:9)
	i1  = [1,2,3,6]
	i2  = [1,4,5,7]
	il1 = [2,3,6]
	il2 = [4,5,7]
	it1 = [1,2,3]
	it2 = [1,4,5]

	# # PD joint controller
	kα = 100.0
	kβ = 10.0
	kp = -[kα, kα, kα, kα, kα, kβ, kβ]
	kd =  0.04*kp
	p.u[1:7] = kp .* (p.q1[iθ] - p.q1[iθ]) + kd .* (qd[iθ])

	# PD joint controller
	kp = -100.0
	kd =  0.04*kp
	p.u[1] = kp * (p.q1[3] - p.qref[3]) + kd * qd[3]


	# 2 feet (heel+toe) in contact
	if all(p.contact)
		kpfx = -200.0
		kdfx = 0.04*kpfx
		kpfz = -400.0
		kdfz = -200.0

		xref = p.q1[1] + p.xdref*p.h_sim
		fx = kpgit stfx*(p.q1[1] - xref) + kdfx*(qd[1] - p.xdref)
		fz = kpfz*(p.q1[2] - p.qref[2]+0.02) + kdfz*qd[2] + model.g*m_flamingo
		f = [fx, fz]
		α = 0.5
		f1 = (1-α)*f
		f2 = α*f
		τ1 = virtual_actuator_torque2(p.model, p.q1, f1, body=:foot_1)
		τ2 = virtual_actuator_torque2(p.model, p.q1, f2, body=:foot_2)
		# τ1[3] = clamp(τ1[3], -20.0, 20.0)
		# τ2[3] = clamp(τ2[3], -20.0, 20.0)
		p.u[il1] .= τ1
		p.u[il2] .= τ2
	end

	# Rescale
	p.u .*= traj.h
    return p.u
end

function detect_contact(model::Flamingo, γ::AbstractVector; contact_threshold=1.5e-2)
	nc = model.dim.c
	contact = [false for i=1:nc]
	for i = 1:nc
		if γ[i] .> contact_threshold
			contact[i] = true
		else
			contact[i] = false
		end
	end
	return contact
end
