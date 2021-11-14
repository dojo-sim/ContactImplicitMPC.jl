"""
    WalledCartpole
	This is cartpole model with two rigid walls constraining the movement of the tip of the cartpole.

	    |                              O mt             |
	    |                               |    |           |
	WALL|                              l | θ |          |WALL
	    |                                 |<-|          |
	    |                                  | |          |
	    ______________________|_____________0|___________
                               <------------>
						               x
	    <----------------------------------------------->
	                          2w
"""
mutable struct WalledCartpole{T} <: Model{T}
    nq::Int 
	nu::Int
	nw::Int
	nc::Int

	mb::T # mass of the base of the cartpole
	mt::T # mass of the tip of the cartpole
	mw::T # mass of the wall
	l::T  # length of the pole
	lc::T # length between base and pole center of mass
	w::T  # width between the two walls is
	k::T  # wall stiffness

    μ_world::T # coefficient of friction
    μ_joint::T #
	g::T # gravity

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end


# Kinematics
function _kinematics(model::WalledCartpole, q; mode = :tip)
	θ, x  = q
	if mode == :tip
		return [x - 1.0 * model.l * sin(θ);
				    1.0 * model.l * cos(θ)]
	elseif mode == :base
		return [x, 0.0]
	else
		@error "incorrect mode"
	 	return zeros(2)
	end
end

function _jacobian(model::WalledCartpole, q; mode = :tip)
	θ, x = q

	if mode == :tip
 		return [-1.0 * model.l * cos(θ) 1.0 0.0 0.0;
		        -1.0 * model.l * sin(θ) 0.0 0.0 0.0]
	elseif mode == :base
		return [0.0 1.0 0.0 0.0;
		        0.0 0.0 0.0 0.0]
 	else
 		@error "incorrect mode"
 	 	return
 	end
end

function kinematics(model::WalledCartpole, q)
	[_kinematics(model, q, mode = :tip);
	 _kinematics(model, q, mode = :tip)]
end

function lagrangian(model::WalledCartpole, q, q̇)
	# https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-832-underactuated-robotics-spring-2009/readings/MIT6_832s09_read_ch03.pdf
	mb = model.mb
	mt = model.mt
	mw = model.mw
	lc = model.lc
	k = model.k

	θ = q[1]
	x = q[2]
	xw1 = q[3]
	xw2 = q[4]

	ω = q̇[1]
	ẋ = q̇[2]
	ẋw1 = q̇[3]
	ẋw2 = q̇[4]

	T = 0.5 * (mt + mb) * ẋ^2 - mt * ẋ * ω * lc * cos(θ) + 0.5 * mt * lc^2 * ω^2
	T += 0.5 * mw * ẋw1^2 + 0.5 * mw * ẋw2^2
	V = mt * g * lc * cos(θ)
	V += k * xw1^2 + k * xw2^2
	L = T - V
	return L
end

function ϕ_func(model::WalledCartpole, env::Environment, q)
	# walls at x = -w + xw1, x = w + xw2
	# + x - xw1 + w > 0
	# - x + xw2 + w > 0
	xw1 = q[3]
	xw2 = q[4]
    SVector{2}([_kinematics(model, q, mode = :tip)[1] - xw1 + model.w;
	            model.w + xw2 - _kinematics(model, q, mode = :tip)[1]])
end

function J_func(model::WalledCartpole, env::Environment, q)
	r1 = [0.0 -1.0; 1.0 0.0]
	r2 = [0.0 1.0; -1.0 0.0]

    SMatrix{4, 4}([r1 * (_jacobian(model, q, mode = :tip) + [0 0 -1.  0.; 0 0 0. 0.]);
		           r2 * (_jacobian(model, q, mode = :tip) + [0 0  0. -1.; 0 0 0. 0.])])
end

function B_func(model::WalledCartpole, q)
	@SMatrix [0.0 1.0 0.0 0.0]
end

function A_func(::WalledCartpole, q)
	@SMatrix [1.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0;
			  0.0 0.0 1.0 0.0;
			  0.0 0.0 0.0 1.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 0.1  # coefficient of friction
μ_joint = 1.0

mb = 0.978 # base mass
mt = 0.411 # mass of tip
mw = 0.1 # mass of the wall
l = 0.6 # leg length
lc = 0.4267  # length between the base and the center of mass of the pole
w = 0.35 # half width between walls
k = 50.0 # wall stiffness

# Dimensions
nq = 4
nu = 1
nw = nq
nc = 2
nquat = 0

walledcartpole = WalledCartpole(nq,nu,nw,nc,
			   mb, mt, mw, l, lc, w, k,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(),
			   SVector{nq}(μ_joint * [0.0; 1.0; 3.0; 3.0]))

function friction_coefficients(model::WalledCartpole) 
	return [model.μ_world]
end