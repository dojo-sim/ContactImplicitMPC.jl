"""
    Quadruped 3D
        orientation representation: modified rodrigues parameters
		s = (px, py, pz,
			 r1, r2, r3,
			 αs1, αt1, αc1,
			 αs2, αt2, αc2,
			 αs3, αt3, αc3,
			 αs4, αt4, αc4)
		  = (3d_position,
		  	 MRP,
			 shoulder1, thigh1, calf1, |-> Absolute angle representation with
			 shoulder2, thigh2, calf2, |-> in the body frame.
			 shoulder3, thigh3, calf3, |-> We cannot do full absolute angle
			 shoulder4, thigh4, calf4, |-> representation to respect the joints.
			 )
		The body frame defined by 3 vectors expressed in the world frame:
			R = MRP(rot) = [xb | yb | zb]_W
"""
mutable struct Quadruped3D{T} <: ContactModel
	dim::Dimensions

	g::T
	μ_world::T
	μ_joint::T

	orientation::Symbol
	shoulder_lateral_offset::T

    # torso
    l_torso::T
    d_torso::T
    m_torso::T
    J_torso::T

    # leg 1
		# shoulder
	l_shoulder1::T
	d_shoulder1::T
	m_shoulder1::T
	J_shoulder1::T
        # thigh
    l_thigh1::T
    d_thigh1::T
    m_thigh1::T
    J_thigh1::T
        # calf
    l_calf1::T
    d_calf1::T
    m_calf1::T
    J_calf1::T

    # leg 2
		# shoulder
	l_shoulder2::T
	d_shoulder2::T
	m_shoulder2::T
	J_shoulder2::T
        # thigh
    l_thigh2::T
    d_thigh2::T
    m_thigh2::T
    J_thigh2::T
        # calf
    l_calf2::T
    d_calf2::T
    m_calf2::T
    J_calf2::T

	# leg 3
		# shoulder
	l_shoulder3::T
	d_shoulder3::T
	m_shoulder3::T
	J_shoulder3::T
        # thigh
    l_thigh3::T
    d_thigh3::T
    m_thigh3::T
    J_thigh3::T
        # calf
    l_calf3::T
    d_calf3::T
    m_calf3::T
    J_calf3::T

	# leg 4
		# shoulder
	l_shoulder4::T
	d_shoulder4::T
	m_shoulder4::T
	J_shoulder4::T
        # thigh
    l_thigh4::T
    d_thigh4::T
    m_thigh4::T
    J_thigh4::T
        # calf
    l_calf4::T
    d_calf4::T
    m_calf4::T
    J_calf4::T

	# fast methods
	base
	dyn
	con
	res
	linearized

	spa::SparseStructure

	joint_friction

	env::Environment
end

# Trunk model
#                             middle   com
#                                267 mm
#             o---------------------|---X-----------------o
#                                   13mm
# BRhip o___________________________|___________________________o FRhip
#                    183 mm                    183mm
# kinematics
function kinematics_1(model::Quadruped3D, q; body = :torso, mode = :ee)
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

function jacobian_1(model::Quadruped3D, q; body = :torso, mode = :ee)
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

function kinematics_2(model::Quadruped3D, q; body = :thigh_1, mode = :ee)
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

function jacobian_2(model::Quadruped3D, q; body = :thigh_1, mode = :ee)

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

function kinematics_3(model::Quadruped3D, q; body = :calf_1, mode = :ee)
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

function jacobian_3(model::Quadruped3D, q; body = :calf_1, mode = :ee)

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

function kinematics_4(model::Quadruped3D, q; body = :calf_3, mode = :ee)
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

function jacobian_4(model::Quadruped3D, q; body = :calf_3, mode = :ee)

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


# Lagrangian
function lagrangian(model::Quadruped3D, q, q̇)
	L = 0.0

	# torso
	p_torso = kinematics_1(model, q, body = :torso, mode = :com)
	J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	v_torso = J_torso * q̇

	L += 0.5 * model.m_torso * transpose(v_torso) * v_torso
	L += 0.5 * model.J_torso * q̇[5]^2.0 # wrong, we need to compute ω, and J needs to be a 3x3 matrix
	L -= model.m_torso * model.g * p_torso[3]

	# shoulder 1
	p_shoulder_1 = kinematics_1(model, q, body = :shoulder_1, mode = :com)
	J_shoulder_1 = jacobian_1(model, q, body = :shoulder_1, mode = :com)
	v_shoulder_1 = J_shoulder_1 * q̇

	L += 0.5 * model.m_shoulder1 * transpose(v_shoulder_1) * v_shoulder_1
	L += 0.5 * model.J_shoulder1 * q̇[7]^2.0
	L -= model.m_shoulder1 * model.g * p_shoulder_1[3]

	# shoulder 2
	p_shoulder_2 = kinematics_1(model, q, body = :shoulder_2, mode = :com)
	J_shoulder_2 = jacobian_1(model, q, body = :shoulder_2, mode = :com)
	v_shoulder_2 = J_shoulder_2 * q̇

	L += 0.5 * model.m_shoulder2 * transpose(v_shoulder_2) * v_shoulder_2
	L += 0.5 * model.J_shoulder2 * q̇[10]^2.0
	L -= model.m_shoulder2 * model.g * p_shoulder_2[3]

	# shoulder 3
	p_shoulder_3 = kinematics_2(model, q, body = :shoulder_3, mode = :com)
	J_shoulder_3 = jacobian_2(model, q, body = :shoulder_3, mode = :com)
	v_shoulder_3 = J_shoulder_3 * q̇

	L += 0.5 * model.m_shoulder3 * transpose(v_shoulder_3) * v_shoulder_3
	L += 0.5 * model.J_shoulder3 * q̇[13]^2.0
	L -= model.m_shoulder3 * model.g * p_shoulder_3[3]

	# shoulder 4
	p_shoulder_4 = kinematics_2(model, q, body = :shoulder_4, mode = :com)
	J_shoulder_4 = jacobian_2(model, q, body = :shoulder_4, mode = :com)
	v_shoulder_4 = J_shoulder_4 * q̇

	L += 0.5 * model.m_shoulder4 * transpose(v_shoulder_4) * v_shoulder_4
	L += 0.5 * model.J_shoulder4 * q̇[16]^2.0
	L -= model.m_shoulder4 * model.g * p_shoulder_4[3]


	# thigh 1
	p_thigh_1 = kinematics_2(model, q, body = :thigh_1, mode = :com)
	J_thigh_1 = jacobian_2(model, q, body = :thigh_1, mode = :com)
	v_thigh_1 = J_thigh_1 * q̇

	L += 0.5 * model.m_thigh1 * transpose(v_thigh_1) * v_thigh_1
	L += 0.5 * model.J_thigh1 * q̇[8]^2.0
	L -= model.m_thigh1 * model.g * p_thigh_1[3]

	# thigh 2
	p_thigh_2 = kinematics_2(model, q, body = :thigh_2, mode = :com)
	J_thigh_2 = jacobian_2(model, q, body = :thigh_2, mode = :com)
	v_thigh_2 = J_thigh_2 * q̇

	L += 0.5 * model.m_thigh2 * transpose(v_thigh_2) * v_thigh_2
	L += 0.5 * model.J_thigh2 * q̇[11]^2.0
	L -= model.m_thigh2 * model.g * p_thigh_2[3]

	# thigh 3
	p_thigh_3 = kinematics_3(model, q, body = :thigh_3, mode = :com)
	J_thigh_3 = jacobian_3(model, q, body = :thigh_3, mode = :com)
	v_thigh_3 = J_thigh_3 * q̇

	L += 0.5 * model.m_thigh3 * transpose(v_thigh_3) * v_thigh_3
	L += 0.5 * model.J_thigh3 * q̇[14]^2.0
	L -= model.m_thigh3 * model.g * p_thigh_3[3]

	# thigh 4
	p_thigh_4 = kinematics_3(model, q, body = :thigh_4, mode = :com)
	J_thigh_4 = jacobian_3(model, q, body = :thigh_4, mode = :com)
	v_thigh_4 = J_thigh_4 * q̇

	L += 0.5 * model.m_thigh4 * transpose(v_thigh_4) * v_thigh_4
	L += 0.5 * model.J_thigh4 * q̇[17]^2.0
	L -= model.m_thigh4 * model.g * p_thigh_4[3]


	# calf 1
	p_calf_1 = kinematics_3(model, q, body = :calf_1, mode = :com)
	J_calf_1 = jacobian_3(model, q, body = :calf_1, mode = :com)
	v_calf_1 = J_calf_1 * q̇

	L += 0.5 * model.m_calf1 * transpose(v_calf_1) * v_calf_1
	L += 0.5 * model.J_calf1 * q̇[8]^2.0
	L -= model.m_calf1 * model.g * p_calf_1[3]

	# calf 2
	p_calf_2 = kinematics_3(model, q, body = :calf_2, mode = :com)
	J_calf_2 = jacobian_3(model, q, body = :calf_2, mode = :com)
	v_calf_2 = J_calf_2 * q̇

	L += 0.5 * model.m_calf2 * transpose(v_calf_2) * v_calf_2
	L += 0.5 * model.J_calf2 * q̇[11]^2.0
	L -= model.m_calf2 * model.g * p_calf_2[3]

	# calf 3
	p_calf_3 = kinematics_4(model, q, body = :calf_3, mode = :com)
	J_calf_3 = jacobian_4(model, q, body = :calf_3, mode = :com)
	v_calf_3 = J_calf_3 * q̇

	L += 0.5 * model.m_calf3 * transpose(v_calf_3) * v_calf_3
	L += 0.5 * model.J_calf3 * q̇[14]^2.0
	L -= model.m_calf3 * model.g * p_calf_3[3]

	# calf 4
	p_calf_4 = kinematics_4(model, q, body = :calf_4, mode = :com)
	J_calf_4 = jacobian_4(model, q, body = :calf_4, mode = :com)
	v_calf_4 = J_calf_4 * q̇

	L += 0.5 * model.m_calf4 * transpose(v_calf_4) * v_calf_4
	L += 0.5 * model.J_calf4 * q̇[17]^2.0
	L -= model.m_calf4 * model.g * p_calf_4[3]

	return L
end

function _dLdq(model::Quadruped3D, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::Quadruped3D, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function kinematics(model::Quadruped3D, q)
	p_calf_1 = kinematics_3(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_3(model, q, body = :calf_2, mode = :ee)
	p_calf_3 = kinematics_4(model, q, body = :calf_3, mode = :ee)
	p_calf_4 = kinematics_4(model, q, body = :calf_4, mode = :ee)

	SVector{12}([p_calf_1; p_calf_2; p_calf_3; p_calf_4])
end

# Methods
function M_func(model::Quadruped3D, q)
	M = Diagonal([0.0, 0.0, 0.0,
		0.0, model.J_torso, 0.0,
		model.J_shoulder1, model.J_thigh1, model.J_calf1,
		model.J_shoulder2, model.J_thigh2, model.J_calf2,
		model.J_shoulder3, model.J_thigh3, model.J_calf3,
		model.J_shoulder4, model.J_thigh4, model.J_calf4])

	# torso
	J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	M += model.m_torso * transpose(J_torso) * J_torso

	# shoulder 1
	J_shoulder_1 = jacobian_1(model, q, body = :shoulder_1, mode = :com)
	M += model.m_shoulder1 * transpose(J_shoulder_1) * J_shoulder_1
	# shoulder 2
	J_shoulder_2 = jacobian_1(model, q, body = :shoulder_2, mode = :com)
	M += model.m_shoulder2 * transpose(J_shoulder_2) * J_shoulder_2
	# shoulder 3
	J_shoulder_3 = jacobian_2(model, q, body = :shoulder_3, mode = :com)
	M += model.m_shoulder3 * transpose(J_shoulder_3) * J_shoulder_3
	# shoulder 4
	J_shoulder_4 = jacobian_2(model, q, body = :shoulder_4, mode = :com)
	M += model.m_shoulder4 * transpose(J_shoulder_4) * J_shoulder_4

	# thigh 1
	J_thigh_1 = jacobian_2(model, q, body = :thigh_1, mode = :com)
	M += model.m_thigh1 * transpose(J_thigh_1) * J_thigh_1
	# thigh 2
	J_thigh_2 = jacobian_2(model, q, body = :thigh_2, mode = :com)
	M += model.m_thigh2 * transpose(J_thigh_2) * J_thigh_2
	# thigh 3
	J_thigh_3 = jacobian_3(model, q, body = :thigh_3, mode = :com)
	M += model.m_thigh3 * transpose(J_thigh_3) * J_thigh_3
	# thigh 4
	J_thigh_4 = jacobian_3(model, q, body = :thigh_4, mode = :com)
	M += model.m_thigh4 * transpose(J_thigh_4) * J_thigh_4

	# calf 1
	J_calf_1 = jacobian_3(model, q, body = :calf_1, mode = :com)
	M += model.m_calf1 * transpose(J_calf_1) * J_calf_1
	# calf 2
	J_calf_2 = jacobian_3(model, q, body = :calf_2, mode = :com)
	M += model.m_calf2 * transpose(J_calf_2) * J_calf_2
	# calf 3
	J_calf_3 = jacobian_4(model, q, body = :calf_3, mode = :com)
	M += model.m_calf3 * transpose(J_calf_3) * J_calf_3
	# calf 4
	J_calf_4 = jacobian_4(model, q, body = :calf_4, mode = :com)
	M += model.m_calf4 * transpose(J_calf_4) * J_calf_4
	return M
end

function ϕ_func(model::Quadruped3D, q)
	p_calf_1 = kinematics_3(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_3(model, q, body = :calf_2, mode = :ee)
	p_calf_3 = kinematics_4(model, q, body = :calf_3, mode = :ee)
	p_calf_4 = kinematics_4(model, q, body = :calf_4, mode = :ee)

	SVector{4}([p_calf_1[3] - model.env.surf(p_calf_1[1:1]);
				p_calf_2[3] - model.env.surf(p_calf_2[1:1]);
				p_calf_3[3] - model.env.surf(p_calf_3[1:1]);
				p_calf_4[3] - model.env.surf(p_calf_4[1:1])])
end

function B_func(model::Quadruped3D, q) #wrong, we need to use the matrix MRP(r)
	@SMatrix [
			  0.0  0.0  0.0  -1.0  0.0  0.0   1.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0;
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
			  0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0  0.0  0.0   0.0 -1.0  1.0;
			  ]
end

function A_func(model::Quadruped3D, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Quadruped3D, q)
	J_calf_1 = jacobian_3(model, q, body = :calf_1, mode = :ee)
	J_calf_2 = jacobian_3(model, q, body = :calf_2, mode = :ee)
	J_calf_3 = jacobian_4(model, q, body = :calf_3, mode = :ee)
	J_calf_4 = jacobian_4(model, q, body = :calf_4, mode = :ee)

	return [J_calf_1;
			J_calf_2;
			J_calf_3;
			J_calf_4]
end

function C_func(model::Quadruped3D, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function contact_forces(model::Quadruped3D, γ1, b1, q2, k)
	k = kinematics(model, q2)
	m = friction_mapping(model.env)

	SVector{16}([transpose(rotation(model.env, k[1:1])) * [m * b1[1:4];   γ1[1]];
				 transpose(rotation(model.env, k[3:3])) * [m * b1[5:8];   γ1[2]];
				 transpose(rotation(model.env, k[5:5])) * [m * b1[9:12];  γ1[3]];
				 transpose(rotation(model.env, k[7:7])) * [m * b1[13:16]; γ1[4]]])
end

function velocity_stack(model::Quadruped3D, q1, q2, k, h)
	k = kinematics(model, q2)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(model.env, k[1:1]) * v[1:3]
	v2_surf = rotation(model.env, k[3:3]) * v[4:6]
	v3_surf = rotation(model.env, k[5:5]) * v[7:9]
	v4_surf = rotation(model.env, k[7:7]) * v[10:12]

	v1_surf = v[1:2]
	v2_surf = v[3:4]
	v3_surf = v[5:6]
	v4_surf = v[7:8]

	SVector{8}([v1_surf[1:2]; -v1_surf[1:2];
				v2_surf[1:2]; -v2_surf[1:2];
				v3_surf[1:2]; -v3_surf[1:2];
				v4_surf[1:2]; -v4_surf[1:2]])
end


################################################################################
# Instantiation
################################################################################
# Dimensions
nq = 3 + 3 + 4*3          # configuration dimension
nu = 4 + 4 + 4            # control dimension
nc = 4                    # number of contact points
nw = 3
nquat = 0

# World parameters
g = 9.81      # gravity
μ_world = 1.0 # coefficient of friction
μ_joint = 0.1 # coefficient of torque friction at the joints

orientation = :MRP
shoulder_lateral_offset = 0.047

# ~Unitree A1
# Model parameters
m_torso = 4.713
m_shoulder = 0.696
m_thigh = 1.013
m_leg = 0.166

J_torso = 0.01683
J_shoulder = 0.000469246
J_thigh = 0.00552
J_leg = 0.00299

l_torso = 0.183*2
l_shoulder = 0.08505
l_thigh = 0.2
l_leg = 0.2

d_torso = 0.5 * l_torso + 0.0127
d_shoulder = 0.5 * l_shoulder -0.003311
d_thigh = 0.5 * l_thigh - 0.00323
d_leg = 0.5 * l_leg - 0.006435

m_payload = 5.0
J_payload = 0.05

quadruped_3D = Quadruped3D(Dimensions(nq, nu, nw, nc, nquat),
				g, μ_world, μ_joint,
				orientation, shoulder_lateral_offset,
				l_torso, d_torso, m_torso, J_torso,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				BaseMethods(), DynamicsMethods(), ContactMethods(),
				ResidualMethods(), ResidualMethods(),
				SparseStructure(spzeros(0, 0), spzeros(0, 0)),
				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]),
				environment_2D_flat())

quadruped_payload_3D = Quadruped3D(Dimensions(nq, nu, nw, nc, nquat),
				g, μ_world, μ_joint,
				orientation, shoulder_lateral_offset,
				l_torso, d_torso,
				m_torso + m_payload,
				J_torso + J_payload,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_shoulder, d_shoulder, m_shoulder, J_shoulder,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				BaseMethods(), DynamicsMethods(), ContactMethods(),
				ResidualMethods(), ResidualMethods(),
				SparseStructure(spzeros(0, 0), spzeros(0, 0)),
				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]),
				environment_2D_flat())
