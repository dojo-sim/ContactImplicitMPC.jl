function mrp_quaternion_map(mrp)
	n2 = mrp[1] * mrp[1] + mrp[2] * mrp[2] + mrp[3] * mrp[3]
    M = 2/(1+n2)
    return SVector{4}([(1-n2)/(1+n2), M*mrp[1], M*mrp[2], M*mrp[3]])
end

mrp_rotation_matrix(mrp) = quaternion_rotation_matrix(mrp_quaternion_map(mrp))

mrp = rand(3)

norm(MRP(mrp...) - mrp_rotation_matrix(mrp))

@variables mrp_sym[1:3]


MRP(mrp_sym...)




mrp_rotation_matrix(mrp_sym)
