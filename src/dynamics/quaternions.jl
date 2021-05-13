function conjugate(q)
	s = q[1]
	v = q[2:4]

	return [s; -v]
end

function L_multiply(q)
	s = q[1]
	v = q[2:4]


	SMatrix{4,4}([s -transpose(v);
	              v s * I + skew(v)])
end

function R_multiply(q)
	s = q[1]
	v = q[2:4]


	SMatrix{4,4}([s -transpose(v);
	              v s * I - skew(v)])
end

function multiply(q1, q2)
	L_multiply(q1) * q2
end

# eq. 14 http://roboticexplorationlab.org/papers/planning_with_attitude.pdf
function attitude_jacobian(q)
	s = q[1]
	v = q[2:4]

	[-transpose(v);
	 s * I + skew(v)]
end

# eq. 24 http://roboticexplorationlab.org/papers/Variational_Integrator.pdf
function ω_finite_difference(q1, q2, h)
	2.0 * transpose(attitude_jacobian(q1)) * multiply(conjugate(q1), (q2 - q1) ./ h)
end

# https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
function quaternion_rotation_matrix(q)
	r, i, j, k  = q

	r11 = 1.0 - 2.0 * (j^2.0 + k^2.0)
	r12 = 2.0 * (i * j - k * r)
	r13 = 2.0 * (i * k + j * r)

	r21 = 2.0 * (i * j + k * r)
	r22 = 1.0 - 2.0 * (i^2.0 + k^2.0)
	r23 = 2.0 * (j * k - i * r)

	r31 = 2.0 * (i * k - j * r)
	r32 = 2.0 * (j * k + i * r)
	r33 = 1.0 - 2.0 * (i^2.0 + j^2.0)

	SMatrix{3,3}([r11 r12 r13;
	              r21 r22 r23;
				  r31 r32 r33])
end
