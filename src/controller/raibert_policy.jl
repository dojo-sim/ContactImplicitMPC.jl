"""
    Raibert hopper policy for Hopper2D
"""

@with_kw mutable struct RaibertOptions
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct Raibert <: Policy
	kr_c::Real
	kr_p::Real
	kr_v_stance::Real
	kr_v_flight::Real
	kθ_c::Real
	kθ_p::Real
	kθ_v::Real
	u::AbstractVector
	contact::Bool
	v0::Real
	Tstance::Real
	Tflight::Real
	q0::AbstractVector
	q1::AbstractVector
	opts::RaibertOptions
end

function raibert_policy(model::Hopper2D; v0=0.5, Tstance=0.13, Tflight=0.62,
		kr_c =  8e1,
		kr_p = -1e3,
		kr_v_stance = -1e-2,
		kr_v_flight = -1e1,
		kθ_c =  0.0,
		kθ_p = -3e1,
		kθ_v = -1e1,
		mpc_opts = RaibertOptions(),
		)

	h = ref_traj.h
	u  = zeros(model.dim.u)
	q0 = zeros(model.dim.q)
	q1 = zeros(model.dim.q)
	contact = false
	Raibert(kr_c, kr_p, kr_v_stance, kr_v_flight, kθ_c, kθ_p, kθ_v,
		u, contact, v0, Tstance, Tflight, copy(q0), copy(q1), mpc_opts)
end

function policy(p::Raibert, x, traj, t)
	# Initialization
	h = traj.h
	p.q0 = copy(x)
	p.q1 = copy(traj.q[t+1])

	# Detect contact
	if any(traj.γ[max(1,t-1)] .> 1.5e-2)
		p.contact = true
	else
		p.contact = false
	end

	# Velocities
	qv = (p.q1 - p.q0)/h
	θv = qv[3]
	rv = qv[4]

	# References
	θ1 = p.q1[3]
	r1 = p.q1[4]
	rref = 0.5
	# td = touchdown
	θtd  = asin(p.v0*Tstance/(2*rref))/2

	# Gains
	kr_c = p.kr_c
	kr_p = p.kr_p
	kr_v_stance = p.kr_v_stance
	kr_v_flight = p.kr_v_flight

	kθ_c = p.kθ_c
	kθ_p = p.kθ_p
	kθ_v = p.kθ_v

	if p.contact
		# regulate around θtd using same gain kθ_p but taking into account
		# that stance time is much lower than flight time
		p.u[1] = kθ_c + kθ_p*(θ1 + θtd) * p.Tflight/p.Tstance
		p.u[2] = kr_c + kr_p*(r1 - rref) + kr_v_stance*rv
	else
		# regulate around θtd
		p.u[1] = kθ_p*(θ1 - θtd)  + kθ_v*θv
		p.u[2] = kr_p*(r1 - rref) + kr_v_flight*rv
	end
	p.u .*= h
    return p.u
end
