function conjugate(q)
	v = q[1:3]
	s = q[4]

	return [-v; s]
end

function multiply(q1, q2)
	v1 = q1[1:3]
	s1 = q1[4]

	v2 = q2[1:3]
	s2 = q2[4]

	[cross(v1, v2) + s1 .* v2 + s2 .* v1;
	 s1 * s2 .- transpose(v1) * v2]
end

# eq. 14 http://roboticexplorationlab.org/papers/planning_with_attitude.pdf
function attitude_jacobian(q)
	[-transpose(q[1:3]);
	Diagonal(q[4] * ones(3)) + skew(q[1:3])]
end

# eq. 24 http://roboticexplorationlab.org/papers/Variational_Integrator.pdf
function ω_finite_difference(q1, q2, h)
	2.0 * transpose(attitude_jacobian(q1)) * multiply(conjugate(q1), (q2 - q1) ./ h)
end

# https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
function quaternion_rotation_matrix(q)
	i, j, k, r = q

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
