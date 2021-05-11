# kinematics
function kinematics_1(model::Quad32, q; body = :torso, mode = :ee)
	x = q[1]
	y = q[2]
	z = q[3]
	p = [x; y; z]
	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]

	if body == :torso
		l = model.l_torso
		d = model.d_torso
		if mode == :ee
			return p + l*zb
		elseif mode == :com
			return p + d*zb
		else
			@error "incorrect mode specification"
		end
	elseif body == :shoulder_1
		off = model.shoulder_lateral_offset
		l = model.l_shoulder1
		d = model.d_shoulder1
		αs = q[7]
		if mode == :ee
			return p + off*yb + l*cos(αs)*yb - l*sin(αs)*xb
		elseif mode == :com
			return p + off*yb + d*cos(αs)*yb - d*sin(αs)*xb
		else
			@error "incorrect mode specification"
		end
	elseif body == :shoulder_2
		off = model.shoulder_lateral_offset
		l = model.l_shoulder2
		d = model.d_shoulder2
		αs = q[10]
		if mode == :ee
			return p - off*yb - l*cos(αs)*yb - l*sin(αs)*xb
		elseif mode == :com
			return p - off*yb - d*cos(αs)*yb - d*sin(αs)*xb
		else
			@error "incorrect mode specification"
		end
	else
		@error "incorrect body specification"
	end
end

function jacobian_1(model::Quad32, q; body = :torso, mode = :ee)
	jac = zeros(eltype(q), 3, model.dim.q)
	jac[1, 1] = 1.0
	jac[2, 2] = 1.0
	jac[3, 3] = 1.0

	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]
	# https://github.com/JuliaGeometry/Rotations.jl/blob/master/src/mrps.jl
	∇r_xb = Rotations.∇rotate(R, [1,0,0.])
	∇r_yb = Rotations.∇rotate(R, [0,1,0.])
	∇r_zb = Rotations.∇rotate(R, [0,0,1.])

	if body == :torso
		e = mode == :ee ? model.l_torso : model.d_torso
		# p + e*zb
		jac[1:3, 4:6] = e*∇r_zb
	elseif body == :shoulder_1
		off = model.shoulder_lateral_offset
		e = mode == :ee ? model.l_shoulder1 : model.d_shoulder1
		αs = q[7]
		# p + (off + e*cos(αs))*yb - e*sin(αs)*xb
		jac[1:3, 4:6] = (off + e*cos(αs))*∇r_yb - e*sin(αs)*∇r_xb
		jac[1:3, 7] = -e*sin(αs)*yb - e*cos(αs)*xb
	elseif body == :shoulder_2
		off = model.shoulder_lateral_offset
		e = mode == :ee ? model.l_shoulder2 : model.d_shoulder2
		αs = q[10]
		# p + (-off - e*cos(αs))*yb - e*sin(αs)*xb
		jac[1:3, 4:6] = (-off - e*cos(αs))*∇r_yb - e*sin(αs)*∇r_xb
		jac[1:3, 10] = +e*sin(αs)*yb - e*cos(αs)*xb
	else
		@error "incorrect body specification"
	end
	return jac
end

function kinematics_2(model::Quad32, q; body = :thigh_1, mode = :ee)
	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]

	if body == :shoulder_3
		p = kinematics_1(model, q, body=:torso, mode=:ee)
		off = model.shoulder_lateral_offset
		l = model.l_shoulder3
		d = model.d_shoulder3
		αs = q[13]
		if mode == :ee
			return p + off*yb + l*cos(αs)*yb - l*sin(αs)*xb
		elseif mode == :com
			return p + off*yb + d*cos(αs)*yb - d*sin(αs)*xb
		else
			@error "incorrect mode specification"
		end
	elseif body == :shoulder_4
		p = kinematics_1(model, q, body=:torso, mode=:ee)
		off = model.shoulder_lateral_offset
		l = model.l_shoulder4
		d = model.d_shoulder4
		αs = q[16]
		if mode == :ee
			return p - off*yb - l*cos(αs)*yb - l*sin(αs)*xb
		elseif mode == :com
			return p - off*yb - d*cos(αs)*yb - d*sin(αs)*xb
		else
			@error "incorrect mode specification"
		end

	# Define the leg plane where the calf and thigh lie (x_leg, z_leg).
	# Define z_leg = zb the unit vector aligned with the torso (body forward direction zb)and equal to zb.
	# Define x_leg = xl the unit vector "aligned" with the body downward
	# direction (xb) and contained in the leg plane: x_leg = β*xb + γ*yb
	elseif body == :thigh_1
		p = kinematics_1(model, q, body=:shoulder_1, mode=:ee)
		l = model.l_thigh1
		d = model.d_thigh1
		αs = q[7]
		αt = q[8]
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αt)*xl + l*sin(αt)*zl
		elseif mode == :com
			return p + d*cos(αt)*xl + d*sin(αt)*zl
		else
			@error "incorrect mode specification"
		end
	elseif body == :thigh_2
		p = kinematics_1(model, q, body=:shoulder_2, mode=:ee)
		l = model.l_thigh2
		d = model.d_thigh2
		αs = q[10]
		αt = q[11]
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αt)*xl + l*sin(αt)*zl
		elseif mode == :com
			return p + d*cos(αt)*xl + d*sin(αt)*zl
		else
			@error "incorrect mode specification"
		end
	else
		@error "incorrect body specification"
	end
end

function jacobian_2(model::Quad32, q; body = :thigh_1, mode = :ee)

	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]
	# https://github.com/JuliaGeometry/Rotations.jl/blob/master/src/mrps.jl
	∇r_xb = Rotations.∇rotate(R, [1,0,0.])
	∇r_yb = Rotations.∇rotate(R, [0,1,0.])
	∇r_zb = Rotations.∇rotate(R, [0,0,1.])


	if body == :shoulder_3
		off = model.shoulder_lateral_offset
		e = mode == :ee ? model.l_shoulder3 : model.d_shoulder3
		αs = q[13]
		# p + (off + e*cos(αs))*yb - e*sin(αs)*xb
		jac = jacobian_1(model, q, body=:torso, mode=:ee)
		jac[1:3, 4:6] += (off + e*cos(αs))*∇r_yb - e*sin(αs)*∇r_xb
		jac[1:3, 13] += -e*sin(αs)*yb - e*cos(αs)*xb
	elseif body == :shoulder_4
		off = model.shoulder_lateral_offset
		e = mode == :ee ? model.l_shoulder4 : model.d_shoulder4
		αs = q[16]
		# p + (-off - e*cos(αs))*yb - e*sin(αs)*xb
		jac = jacobian_1(model, q, body=:torso, mode=:ee)
		jac[1:3, 4:6] += (-off - e*cos(αs))*∇r_yb - e*sin(αs)*∇r_xb
		jac[1:3, 16] += +e*sin(αs)*yb - e*cos(αs)*xb
	elseif body == :thigh_1
		e = mode == :ee ? model.l_thigh1 : model.d_thigh1
		αs = q[7]
		αt = q[8]
		# p + e*cos(αt)*xl + e*sin(αt)*zl
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb + sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_1(model, q, body=:shoulder_1, mode=:ee)
		jac[1:3, 4:6] += e*cos(αt)*∇r_xl + e*sin(αt)*∇r_zl
		jac[1:3, 7] += e*cos(αt) * (-sin(αs)*xb + cos(αs)*yb)
		jac[1:3, 8] += -e*sin(αt)*xl + e*cos(αt)*zl
	elseif body == :thigh_2
		e = mode == :ee ? model.l_thigh2 : model.d_thigh2
		αs = q[10]
		αt = q[11]
		# p + e*cos(αt)*xl + e*sin(αt)*zl
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb - sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_1(model, q, body=:shoulder_2, mode=:ee)
		jac[1:3, 4:6] += e*cos(αt)*∇r_xl + e*sin(αt)*∇r_zl
		jac[1:3, 10] += e*cos(αt) * (-sin(αs)*xb - cos(αs)*yb)
		jac[1:3, 11] += -e*sin(αt)*xl + e*cos(αt)*zl
	else
		@error "incorrect body specification"
	end
	return jac
end

function kinematics_3(model::Quad32, q; body = :calf_1, mode = :ee)
	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]

	if body == :thigh_3
		p = kinematics_2(model, q, body=:shoulder_3, mode=:ee)
		l = model.l_thigh3
		d = model.d_thigh3
		αs = q[13]
		αt = q[14]
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αt)*xl + l*sin(αt)*zl
		elseif mode == :com
			return p + d*cos(αt)*xl + d*sin(αt)*zl
		else
			@error "incorrect mode specification"
		end
	elseif body == :thigh_4
		p = kinematics_2(model, q, body=:shoulder_4, mode=:ee)
		l = model.l_thigh4
		d = model.d_thigh4
		αs = q[16]
		αt = q[17]
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αt)*xl + l*sin(αt)*zl
		elseif mode == :com
			return p + d*cos(αt)*xl + d*sin(αt)*zl
		else
			@error "incorrect mode specification"
		end
	elseif body == :calf_1
		p = kinematics_2(model, q, body=:thigh_1, mode=:ee)
		l = model.l_calf1
		d = model.d_calf1
		αs = q[7]
		αc = q[9]
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αc)*xl + l*sin(αc)*zl
		elseif mode == :com
			return p + d*cos(αc)*xl + d*sin(αc)*zl
		else
			@error "incorrect mode specification"
		end
	elseif body == :calf_2
		p = kinematics_2(model, q, body=:thigh_2, mode=:ee)
		l = model.l_calf2
		d = model.d_calf2
		αs = q[10]
		αc = q[12]
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αc)*xl + l*sin(αc)*zl
		elseif mode == :com
			return p + d*cos(αc)*xl + d*sin(αc)*zl
		else
			@error "incorrect mode specification"
		end
	else
		@error "incorrect body specification"
	end
end

function jacobian_3(model::Quad32, q; body = :calf_1, mode = :ee)

	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]
	# https://github.com/JuliaGeometry/Rotations.jl/blob/master/src/mrps.jl
	∇r_xb = Rotations.∇rotate(R, [1,0,0.])
	∇r_yb = Rotations.∇rotate(R, [0,1,0.])
	∇r_zb = Rotations.∇rotate(R, [0,0,1.])

	if body == :thigh_3
		e = mode == :ee ? model.l_thigh3 : model.d_thigh3
		αs = q[13]
		αt = q[14]
		# p + e*cos(αt)*xl + e*sin(αt)*zl
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb + sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_2(model, q, body=:shoulder_3, mode=:ee)
		jac[1:3, 4:6] += e*cos(αt)*∇r_xl + e*sin(αt)*∇r_zl
		jac[1:3, 13] += e*cos(αt) * (-sin(αs)*xb + cos(αs)*yb)
		jac[1:3, 14] += -e*sin(αt)*xl + e*cos(αt)*zl
	elseif body == :thigh_4
		e = mode == :ee ? model.l_thigh4 : model.d_thigh4
		αs = q[16]
		αt = q[17]
		# p + e*cos(αt)*xl + e*sin(αt)*zl
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb - sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_2(model, q, body=:shoulder_4, mode=:ee)
		jac[1:3, 4:6] += e*cos(αt)*∇r_xl + e*sin(αt)*∇r_zl
		jac[1:3, 16] += e*cos(αt) * (-sin(αs)*xb - cos(αs)*yb)
		jac[1:3, 17] += -e*sin(αt)*xl + e*cos(αt)*zl
	elseif body == :calf_1
		e = mode == :ee ? model.l_calf1 : model.d_calf1
		αs = q[7]
		αc = q[9]
		# p + e*cos(αc)*xl + e*sin(αc)*zl
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb + sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_2(model, q, body=:thigh_1, mode=:ee)
		jac[1:3, 4:6] += e*cos(αc)*∇r_xl + e*sin(αc)*∇r_zl
		jac[1:3, 7] += e*cos(αc) * (-sin(αs)*xb + cos(αs)*yb)
		jac[1:3, 9] += -e*sin(αc)*xl + e*cos(αc)*zl
	elseif body == :calf_2
		e = mode == :ee ? model.l_calf2 : model.d_calf2
		αs = q[10]
		αc = q[12]
		# p + e*cos(αc)*xl + e*sin(αc)*zl
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb - sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_2(model, q, body=:thigh_2, mode=:ee)
		jac[1:3, 4:6] += e*cos(αc)*∇r_xl + e*sin(αc)*∇r_zl
		jac[1:3, 10] += e*cos(αc) * (-sin(αs)*xb - cos(αs)*yb)
		jac[1:3, 12] += -e*sin(αc)*xl + e*cos(αc)*zl
	else
		@error "incorrect body specification"
	end
	return jac
end

function kinematics_4(model::Quad32, q; body = :calf_3, mode = :ee)
	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]

	if body == :calf_3
		p = kinematics_3(model, q, body=:thigh_3, mode=:ee)
		l = model.l_calf3
		d = model.d_calf3
		αs = q[13]
		αc = q[15]
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αc)*xl + l*sin(αc)*zl
		elseif mode == :com
			return p + d*cos(αc)*xl + d*sin(αc)*zl
		else
			@error "incorrect mode specification"
		end
	elseif body == :calf_4
		p = kinematics_3(model, q, body=:thigh_4, mode=:ee)
		l = model.l_calf4
		d = model.d_calf4
		αs = q[16]
		αc = q[18]
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		if mode == :ee
			return p + l*cos(αc)*xl + l*sin(αc)*zl
		elseif mode == :com
			return p + d*cos(αc)*xl + d*sin(αc)*zl
		else
			@error "incorrect mode specification"
		end
	else
		@error "incorrect body specification"
	end
end

function jacobian_4(model::Quad32, q; body = :calf_3, mode = :ee)

	R = eval(model.orientation)(view(q, 4:6)...)
	xb = R[1:3,1]
	yb = R[1:3,2]
	zb = R[1:3,3]
	# https://github.com/JuliaGeometry/Rotations.jl/blob/master/src/mrps.jl
	∇r_xb = Rotations.∇rotate(R, [1,0,0.])
	∇r_yb = Rotations.∇rotate(R, [0,1,0.])
	∇r_zb = Rotations.∇rotate(R, [0,0,1.])

	if body == :calf_3
		e = mode == :ee ? model.l_calf3 : model.d_calf3
		αs = q[13]
		αc = q[15]
		# p + e*cos(αc)*xl + e*sin(αc)*zl
		xl = cos(αs)*xb + sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb + sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_3(model, q, body=:thigh_3, mode=:ee)
		jac[1:3, 4:6] += e*cos(αc)*∇r_xl + e*sin(αc)*∇r_zl
		jac[1:3, 13] += e*cos(αc) * (-sin(αs)*xb + cos(αs)*yb)
		jac[1:3, 15] += -e*sin(αc)*xl + e*cos(αc)*zl
	elseif body == :calf_4
		e = mode == :ee ? model.l_calf4 : model.d_calf4
		αs = q[16]
		αc = q[18]
		# p + e*cos(αc)*xl + e*sin(αc)*zl
		xl = cos(αs)*xb - sin(αs)*yb
		zl = zb
		∇r_xl = cos(αs)*∇r_xb - sin(αs)*∇r_yb
		∇r_zl = ∇r_zb
		jac = jacobian_3(model, q, body=:thigh_4, mode=:ee)
		jac[1:3, 4:6] += e*cos(αc)*∇r_xl + e*sin(αc)*∇r_zl
		jac[1:3, 16] += e*cos(αc) * (-sin(αs)*xb - cos(αs)*yb)
		jac[1:3, 18] += -e*sin(αc)*xl + e*cos(αc)*zl
	else
		@error "incorrect body specification"
	end
	return jac
end

function B_func(model::Quad32, q) #wrong, we need to use the matrix MRP(r)
	@SMatrix [0.0  0.0  0.0  -1.0  0.0  0.0   1.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0;
			  0.0  0.0  0.0   0.0 -1.0  0.0   0.0  1.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0;
			  0.0  0.0  0.0   0.0  0.0  0.0   0.0 -1.0  1.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0;

			  0.0  0.0  0.0   1.0  0.0  0.0   0.0  0.0  0.0   1.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0;
			  0.0  0.0  0.0   0.0 -1.0  0.0   0.0  0.0  0.0   0.0  1.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0;
			  0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0 -1.0  1.0   0.0  0.0  0.0   0.0  0.0  0.0;

			  0.0  0.0  0.0  -1.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   1.0  0.0  0.0   0.0  0.0  0.0;
			  0.0  0.0  0.0   0.0 -1.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  1.0  0.0   0.0  0.0  0.0;
			  0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0 -1.0  1.0   0.0  0.0  0.0;

			  0.0  0.0  0.0   1.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   1.0  0.0  0.0;
			  0.0  0.0  0.0   0.0 -1.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  1.0  0.0;
			  0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0 -1.0  1.0]
end

function A_func(model::Quad32, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Quad32, q)
	J_calf_1 = jacobian_3(model, q, body = :calf_1, mode = :ee)
	J_calf_2 = jacobian_3(model, q, body = :calf_2, mode = :ee)
	J_calf_3 = jacobian_4(model, q, body = :calf_3, mode = :ee)
	J_calf_4 = jacobian_4(model, q, body = :calf_4, mode = :ee)

	return [J_calf_1;
			J_calf_2;
			J_calf_3;
			J_calf_4]
end

function kinematics(model::Quad32, q)
	p_calf_1 = kinematics_3(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_3(model, q, body = :calf_2, mode = :ee)
	p_calf_3 = kinematics_4(model, q, body = :calf_3, mode = :ee)
	p_calf_4 = kinematics_4(model, q, body = :calf_4, mode = :ee)

	SVector{12}([p_calf_1; p_calf_2; p_calf_3; p_calf_4])
end

function ϕ_func(model::Quad32, q)
	p_calf_1 = kinematics_3(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_3(model, q, body = :calf_2, mode = :ee)
	p_calf_3 = kinematics_4(model, q, body = :calf_3, mode = :ee)
	p_calf_4 = kinematics_4(model, q, body = :calf_4, mode = :ee)

	SVector{4}([p_calf_1[3] - model.env.surf(p_calf_1[1:1]);
				p_calf_2[3] - model.env.surf(p_calf_2[1:1]);
				p_calf_3[3] - model.env.surf(p_calf_3[1:1]);
				p_calf_4[3] - model.env.surf(p_calf_4[1:1])])
end
