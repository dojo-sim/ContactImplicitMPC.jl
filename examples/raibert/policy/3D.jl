"""
    Raibert hopper policy for Hopper2D
"""

ContactImplicitMPC.@with_kw mutable struct RaibertOptions
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct Raibert3D <: ContactImplicitMPC.Policy
	kr_c::Real
	kr_p::Real
	kr_v_stance::Real
	kr_v_flight::Real
	kθ_c::Real
	kθ_p::Real
	kθ_v::Real
	u::AbstractVector
	contact::Bool
	v0::Vector{Real}
	Tstance::Real
	Tflight::Real
	q0::AbstractVector
	q1::AbstractVector
	opts::RaibertOptions
end

function raibert_policy(model::ContactImplicitMPC.Hopper3D;
		v0=[0.0; 0.0], Tstance=0.13, Tflight=0.62, h = 0.1,
		kr_c =  8e1,
		kr_p = -1e3,
		kr_v_stance = -1e-2,
		kr_v_flight = -1e1,
		kθ_c =  0.0,
		kθ_p = -6e1,
		kθ_v = -1e1,
		mpc_opts = RaibertOptions(),
		)

	@assert v0[1] == 0.0 || v0[2] == 0.0

	u  = zeros(model.dim.u)
	q0 = zeros(model.dim.q)
	q1 = zeros(model.dim.q)
	contact = false

	Raibert3D(kr_c, kr_p, kr_v_stance, kr_v_flight, kθ_c, kθ_p, kθ_v,
		u, contact, v0, Tstance, Tflight, copy(q0), copy(q1), mpc_opts)
end

function ContactImplicitMPC.policy(p::Raibert3D, x, traj, t)
	# Initialization
	h = traj.h
	p.q0 = copy(traj.q[t])
	p.q1 = copy(x)

	# get Euler angles
	rot0 = RotXYZ(MRP(p.q0[4:6]...))
	eul0 = [rot0.theta1; rot0.theta2; rot0.theta3]

	rot1 = RotXYZ(MRP(p.q1[4:6]...))
	eul1 = [rot1.theta1; rot1.theta2; rot1.theta3]

	# Detect contact
	if any(traj.γ[max(1,t-1)] .> 1.5e-2)
		p.contact = true
	else
		p.contact = false
	end


	# Velocities
	qv = (p.q1 - p.q0)/h
	θv = (eul1 - eul0)/h
	rv = qv[7]

	# References
	θ1 = eul1
	r1 = p.q1[7]
	rref = 0.5

	# td = touchdown
	θtd  = [asin(p.v0[1]*Tstance/(2*rref))/2; asin(p.v0[2]*Tstance/(2*rref))/2]

	# Gains
	kr_c = p.kr_c
	kr_p = p.kr_p
	kr_v_stance = p.kr_v_stance
	kr_v_flight = p.kr_v_flight

	kθ_c = p.kθ_c
	kθ_p = p.kθ_p
	kθ_v = p.kθ_v

	dir = qv[1:2] ./ norm(qv[1:2]) # TODO: change to binary switch

	if p.contact
		# regulate around θtd using same gain kθ_p but taking into account
		# that stance time is much lower than flight time
		p.u[1] = (dir[2] * kθ_c + kθ_p*(θ1[1] + dir[2] * θtd[1]) * p.Tflight/p.Tstance)
		p.u[2] = (dir[1] * kθ_c + kθ_p*(θ1[2] + dir[1] * θtd[2]) * p.Tflight/p.Tstance)
		p.u[3] = kr_c + kr_p*(r1 - rref) + kr_v_stance*rv
	else
		# regulate around θtd
		p.u[1] = (kθ_p*(θ1[1] - dir[2] * θtd[1]))  + kθ_v*θv[1]
		p.u[2] = (kθ_p*(θ1[2] - dir[1] * θtd[2]))  + kθ_v*θv[2]
		p.u[3] = kr_p*(r1 - rref) + kr_v_flight*rv
	end

	p.u .*= h

    return p.u
end
