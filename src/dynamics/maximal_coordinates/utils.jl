function G_func(q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

function static_quaternion(q::UnitQuaternion)
	SVector{4}([q.w, q.x, q.y, q.z])
end

function quaternion_offset(q_parent, q_child)
	L_multiply(conjugate(q_parent)) * q_child
end

function get_quaternion(q)
	return view(q, 4:7)
end

friction_map = SMatrix{2, 4}([1.0 0.0 -1.0 0.0;
                   			  0.0 1.0 0.0 -1.0])
