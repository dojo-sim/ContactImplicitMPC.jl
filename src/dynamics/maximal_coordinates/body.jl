struct Body{T,I}
	mass::T
	inertia::I
	gravity::SVector{3,T}
    id::Symbol
end

function Body(mass::T, inertia::I;
    id = :body, gravity = 9.81) where {T,I}
    return Body(
        mass,
        inertia,
        SVector{3,T}([0.0, 0.0, gravity]),
		id)
end

body1 = Body(1.0, Diagonal(@SVector ones(3)), id = :body1)
body2 = Body(1.0, Diagonal(@SMatrix ones(3, 3)), id = :body2)

function dynamics(mass, inertia, gravity, h, q0, q1, q2, f)

	# unpack
	p0 = q0[1:3]
	quat0 = q0[4:7]

	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	# evalutate at midpoint
	vm1 = (p1 - p0) / h[1]
    vm2 = (p2 - p1) / h[1]
	ω1 = ω_finite_difference(quat0, quat1, h)
	ω2 = ω_finite_difference(quat1, quat2, h)

	# linear
	d_linear = mass * (vm1 - vm2) - h[1] * mass * gravity

	# angular
	d_angular = -1.0 * (inertia * ω2 * sqrt(4.0 / h[1]^2.0 - transpose(ω2) * ω2)
		+ cross(ω2, inertia * ω2)
		- inertia * ω1 * sqrt(4.0 / h[1]^2.0 - transpose(ω1) * ω1)
		+ cross(ω1, inertia * ω1))

	# NOTE: torques require 2x (see 'Linear-Time Variational Integrators in Maximal Coordinates')
	return SVector{6}([d_linear; d_angular]) + SVector{6}([1.0; 1.0; 1.0; 2.0; 2.0; 2.0]) .* f
end

using Symbolics, BenchmarkTools

@variables mass, inertia[1:3, 1:3], gravity[1:3], h[1]
@variables q0[1:7], q1[1:7], q2[1:7]
@variables f1[1:6]

d = dynamics(mass, inertia, gravity, h, q0, q1, q2, f1)

d_func = eval(Symbolics.build_function(d, mass, inertia, gravity, h, q0, q1, q2, f1)[1])
d_func! = eval(Symbolics.build_function(d, mass, inertia, gravity, h, q0, q1, q2, f1)[2])

_mass = 1.0
_inertia = SMatrix{3, 3}(rand(3, 3))
_gravity = SVector{3}([0.0, 0.0, 9.81])
_h = 0.1
_q0 = rand(7)
_q1 = rand(7)
_q2 = rand(7)
_f1 = rand(6)
_d = zeros(6)

@benchmark $_d .= d_func($_mass, $_inertia, $_gravity, $_h, $_q0, $_q1, $_q2, $_f1)
@benchmark d_func!($_d, $_mass, $_inertia, $_gravity, $_h, $_q0, $_q1, $_q2, $_f1)
