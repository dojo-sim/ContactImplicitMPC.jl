struct Body
	mass
	inertia
	kinematics
	gravity
end

friction_map = SMatrix{2, 4}([1.0 0.0 -1.0 0.0;
                   			  0.0 1.0 0.0 -1.0])

function rotation_axes(q_parent, q_child, q_offset)
	# q_parent^-1 q_child q_offset^-1
	(R_multiply(conjugate(q_offset)) * L_multiply(conjugate(q_parent)) * q_child)[2:4]
end

function rotation_angular_velocity(q1_parent, q1_child, q2_parent, q2_child, q_offset, h)
	# q_parent^-1 q_child q_offset^-1
	q1 = (R_multiply(conjugate(q_offset)) * L_multiply(conjugate(q1_parent)) * q1_child)
	q2 = (R_multiply(conjugate(q_offset)) * L_multiply(conjugate(q2_parent)) * q2_child)

	ω = ω_finite_difference(q1, q2, h)

	return ω
end

function get_quaternion(q)
	return q[4:7]
end

function G_func(q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

# floating base
fb_length = 1.0
fb_width = 0.5
fb_height = 0.1
floating_base = Body(Diagonal(1.0 * ones(3)),
					 Diagonal([0.1, 0.1, 0.1]),
					 [[0.5 * fb_length, 0.5 * fb_width, 0.0],
					  [0.5 * fb_length, -0.5 * fb_width, 0.0],
					  [-0.5 * fb_length, 0.5 * fb_width, 0.0],
					  [-0.5 * fb_length, -0.5 * fb_width, 0.0]],
					 [0.0; 0.0; 1.0 * 9.81])

# link
link_length = 0.25
link = Body(Diagonal(0.1 * ones(3)),
			Diagonal([0.01, 0.01, 0.01]),
			[[0.0, 0.0, 0.5 * link_length],
			 [0.0, 0.0, -0.5 * link_length]],
			[0.0; 0.0; 1.0 * 9.81])

nq = 7

@variables q0[1:nq], q1[1:nq], q2[1:nq], r[1:3]

function kinematics(r, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	k1 = p + quaternion_rotation_matrix(quat) * r

	return k1
end

k = kinematics(r, q0)
k_func = eval(Symbolics.build_function(k, r, q0)[1])
∇k = Symbolics.jacobian(k, q0) * G_func(q0)
∇k_func = eval(Symbolics.build_function(∇k, r, q0)[1])

@variables quat_parent[1:4], quat_child[1:4], quat_offset[1:4]
ra = rotation_axes(quat_parent, quat_child, quat_offset)
dra_parent = Symbolics.jacobian(ra, quat_parent) * attitude_jacobian(quat_parent)
dra_child = Symbolics.jacobian(ra, quat_child) * attitude_jacobian(quat_child)

ra_func = eval(Symbolics.build_function(ra, quat_parent, quat_child, quat_offset)[1])
dra_parent_func = eval(Symbolics.build_function(dra_parent, quat_parent, quat_child, quat_offset)[1])
dra_child_func = eval(Symbolics.build_function(dra_child, quat_parent, quat_child, quat_offset)[1])

function dynamics(model::Body, h, q0, q1, q2, f, τ)

	p0 = q0[1:3]
	quat0 = q0[4:7]

	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	# evalutate at midpoint
	ω1 = ω_finite_difference(quat0, quat1, h)
	ω2 = ω_finite_difference(quat1, quat2, h)

	qm1 = [0.5 * (p0 + p1); zeros(4)]
    vm1 = [(p1 - p0) / h[1]; ω1]
    qm2 = [0.5 * (p1 + p2); zeros(4)]
    vm2 = [(p2 - p1) / h[1]; ω2]

	d_linear = model.mass * (vm1[1:3] - vm2[1:3]) - h[1] * model.mass * model.gravity
	d_angular = -1.0 * (model.inertia * ω2 * sqrt(4.0 / h[1]^2.0 - transpose(ω2) * ω2)
		+ cross(ω2, model.inertia * ω2)
		- model.inertia * ω1 * sqrt(4.0 / h[1]^2.0 - transpose(ω1) * ω1)
		+ cross(ω1, model.inertia * ω1))

	d_angular .-= 2.0 * τ
	return [d_linear; d_angular] + f
end

@variables f1[1:6], τ1[1:3], h[1:1]

d_link = dynamics(link, h, q0, q1, q2, f1, τ1)
∇d_link = Symbolics.jacobian(d_link, q2) * G_func(q2)

d_link_func = eval(Symbolics.build_function(d_link, h, q0, q1, q2, f1, τ1)[1])
∇d_link_func = eval(Symbolics.build_function(∇d_link, h, q0, q1, q2, f1, τ1)[1])

d_link_func([1.0], ones(7), ones(7), ones(7), ones(3), ones(3))

d_fb = dynamics(floating_base, h, q0, q1, q2, f1, τ1)
∇d_fb = Symbolics.jacobian(d_fb, q2) * G_func(q2)

d_fb_func = eval(Symbolics.build_function(d_fb, h, q0, q1, q2, f1, τ1)[1])
∇d_fb_func = eval(Symbolics.build_function(∇d_fb, h, q0, q1, q2, f1, τ1)[1])

d_fb_func([1.0], ones(7), ones(7), ones(7), ones(3), ones(3))

links = [link, link, link,
		 link, link, link,
		 link, link, link]
		 # link, link, link]

N = length(links)
n_legs = 3
nf = 3 * n_legs
nr = 3 * n_legs

nc = n_legs # number of contact points
nb = 4 # double parameterized friction
n_contact = nc + nc + nc * (nb + nb) + nc + nc

E = SMatrix{nc, nc * nb}(kron(Diagonal(ones(nc)), ones(1, Int(nb))))
μ_world = 0.5 * ones(nc)

x_axis_mask = [0.0 0.0;
			   1.0 0.0;
			   0.0 1.0]

x_axis_input = [1.0 0.0 0.0;
			    0.0 0.0 0.0;
			    0.0 0.0 0.0]

y_axis_mask = [1.0 0.0;
			   0.0 0.0;
			   0.0 1.0]

y_axis_input = [0.0 0.0 0.0;
			    0.0 1.0 0.0;
			    0.0 0.0 0.0]

z_axis_mask = [1.0 0.0;
			   0.0 1.0;
			   0.0 0.0]

z_axis_input = [0.0 0.0 0.0;
			    0.0 0.0 0.0;
			    0.0 0.0 1.0]

# floating base
z_base = 2 * link_length + 0.0
p0_fb = [0.0; 0.0; z_base]
_quat0_fb = one(UnitQuaternion)
quat0_fb = [Rotations.scalar(_quat0_fb); Rotations.vector(_quat0_fb)...]
q0_fb = [p0_fb; quat0_fb]

p1_fb = [0.0; 0.0; z_base]
quat1_fb = copy(quat0_fb)
q1_fb = [p1_fb; quat1_fb]

θ_leg_shoulder = 0.5 * π + 0.0 * π
θ_leg_elbow = 0.25 * π

## LEG 1
# link 1
p0_l1 = p0_fb + floating_base.kinematics[1] + [0.0; 0.5 * link_length * sin(θ_leg_shoulder); 0.5 * link_length * cos(θ_leg_shoulder)]
_quat0_l1 = UnitQuaternion(AngleAxis(-θ_leg_shoulder, 1.0, 0.0, 0.0))
quat0_l1 = [Rotations.scalar(_quat0_l1); Rotations.vector(_quat0_l1)...]
q0_l1 = [p0_l1; quat0_l1]

p1_l1 = copy(p0_l1)
quat1_l1 = copy(quat0_l1)
q1_l1 = [p1_l1; quat1_l1]

rot_off_fb_l1 = L_multiply(conjugate(quat1_fb)) * quat1_l1

# link 2
p0_l2 = k_func(links[1].kinematics[1], q1_l1) + [-0.5 * link_length * sin(θ_leg_elbow); 0.0; -0.5 * link_length * cos(θ_leg_elbow)]
_quat0_l2 = UnitQuaternion(AngleAxis(θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l2 = [Rotations.scalar(_quat0_l2); Rotations.vector(_quat0_l2)...]
q0_l2 = [p0_l2; quat0_l2]

p1_l2 = copy(p0_l2)
quat1_l2 = copy(quat0_l2)
q1_l2 = [p1_l2; quat1_l2]

rot_off_l1_l2 = L_multiply(conjugate(quat1_l1)) * quat1_l2

# link 3
p0_l3 = p0_l2 + [0.0; 0.0; -1.0 * link_length * cos(θ_leg_elbow)]
_quat0_l3 = UnitQuaternion(AngleAxis(-θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l3 = [Rotations.scalar(_quat0_l3); Rotations.vector(_quat0_l3)...]
q0_l3 = [p0_l3; quat0_l3]

p1_l3 = copy(p0_l3)
quat1_l3 = copy(quat0_l3)
q1_l3 = [p1_l3; quat1_l3]

rot_off_l2_l3 = L_multiply(conjugate(quat1_l2)) * quat1_l3

## LEG 2
# link 4
p0_l4 = p0_fb + floating_base.kinematics[2] + [0.0; -0.5 * link_length * sin(θ_leg_shoulder); 0.5 * link_length * cos(θ_leg_shoulder)]
_quat0_l4 = UnitQuaternion(AngleAxis(θ_leg_shoulder, 1.0, 0.0, 0.0))
quat0_l4 = [Rotations.scalar(_quat0_l4); Rotations.vector(_quat0_l4)...]
q0_l4 = [p0_l4; quat0_l4]

p1_l4 = copy(p0_l4)
quat1_l4 = copy(quat0_l4)
q1_l4 = [p1_l4; quat1_l4]

rot_off_fb_l4 = L_multiply(conjugate(quat1_fb)) * quat1_l4

# link 5
p0_l5 = k_func(links[1].kinematics[1], q1_l4) + [-0.5 * link_length * sin(θ_leg_elbow); 0.0; -0.5 * link_length * cos(θ_leg_elbow)]
_quat0_l5 = UnitQuaternion(AngleAxis(θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l5 = [Rotations.scalar(_quat0_l5); Rotations.vector(_quat0_l5)...]
q0_l5 = [p0_l5; quat0_l5]

p1_l5 = copy(p0_l5)
quat1_l5 = copy(quat0_l5)
q1_l5 = [p1_l5; quat1_l5]

rot_off_l4_l5 = L_multiply(conjugate(quat1_l4)) * quat1_l5

# link 6
p0_l6 = p0_l5 + [0.0; 0.0; -1.0 * link_length * cos(θ_leg_elbow)]
_quat0_l6 = UnitQuaternion(AngleAxis(-θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l6 = [Rotations.scalar(_quat0_l6); Rotations.vector(_quat0_l6)...]
q0_l6 = [p0_l6; quat0_l6]

p1_l6 = copy(p0_l6)
quat1_l6 = copy(quat0_l6)
q1_l6 = [p1_l6; quat1_l6]

rot_off_l5_l6 = L_multiply(conjugate(quat1_l5)) * quat1_l6

## LEG 3
# link 7
p0_l7 = p0_fb + floating_base.kinematics[3] + [0.0; 0.5 * link_length * sin(θ_leg_shoulder); 0.5 * link_length * cos(θ_leg_shoulder)]
_quat0_l7 = UnitQuaternion(AngleAxis(-θ_leg_shoulder, 1.0, 0.0, 0.0))
quat0_l7 = [Rotations.scalar(_quat0_l7); Rotations.vector(_quat0_l7)...]
q0_l7 = [p0_l7; quat0_l7]

p1_l7 = copy(p0_l7)
quat1_l7 = copy(quat0_l7)
q1_l7 = [p1_l7; quat1_l7]

rot_off_fb_l7 = L_multiply(conjugate(quat1_fb)) * quat1_l7

# link 8
p0_l8 = k_func(links[1].kinematics[1], q1_l7) + [-0.5 * link_length * sin(θ_leg_elbow); 0.0; -0.5 * link_length * cos(θ_leg_elbow)]
_quat0_l8 = UnitQuaternion(AngleAxis(θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l8 = [Rotations.scalar(_quat0_l8); Rotations.vector(_quat0_l8)...]
q0_l8 = [p0_l8; quat0_l8]

p1_l8 = copy(p0_l8)
quat1_l8 = copy(quat0_l8)
q1_l8 = [p1_l8; quat1_l8]

rot_off_l7_l8 = L_multiply(conjugate(quat1_l7)) * quat1_l8

# link 9
p0_l9 = p0_l8 + [0.0; 0.0; -1.0 * link_length * cos(θ_leg_elbow)]
_quat0_l9 = UnitQuaternion(AngleAxis(-θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l9 = [Rotations.scalar(_quat0_l9); Rotations.vector(_quat0_l9)...]
q0_l9 = [p0_l9; quat0_l9]

p1_l9 = copy(p0_l9)
quat1_l9 = copy(quat0_l9)
q1_l9 = [p1_l9; quat1_l9]

rot_off_l8_l9 = L_multiply(conjugate(quat1_l8)) * quat1_l9

## LEG 4
# link 10
p0_l10 = p0_fb + floating_base.kinematics[4] + [0.0; -0.5 * link_length * sin(θ_leg_shoulder); 0.5 * link_length * cos(θ_leg_shoulder)]
_quat0_l10 = UnitQuaternion(AngleAxis(θ_leg_shoulder, 1.0, 0.0, 0.0))
quat0_l10 = [Rotations.scalar(_quat0_l10); Rotations.vector(_quat0_l10)...]
q0_l10 = [p0_l10; quat0_l10]

p1_l10 = copy(p0_l10)
quat1_l10 = copy(quat0_l10)
q1_l10 = [p1_l10; quat1_l10]

rot_off_fb_l10 = L_multiply(conjugate(quat1_fb)) * quat1_l10

# link 11
p0_l11 = k_func(links[1].kinematics[1], q1_l10) + [-0.5 * link_length * sin(θ_leg_elbow); 0.0; -0.5 * link_length * cos(θ_leg_elbow)]
_quat0_l11 = UnitQuaternion(AngleAxis(θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l11 = [Rotations.scalar(_quat0_l11); Rotations.vector(_quat0_l11)...]
q0_l11 = [p0_l11; quat0_l11]

p1_l11 = copy(p0_l11)
quat1_l11 = copy(quat0_l11)
q1_l11 = [p1_l11; quat1_l11]

rot_off_l10_l11 = L_multiply(conjugate(quat1_l10)) * quat1_l11

# link 12
p0_l12 = p0_l11 + [0.0; 0.0; -1.0 * link_length * cos(θ_leg_elbow)]
_quat0_l12 = UnitQuaternion(AngleAxis(-θ_leg_elbow, 0.0, 1.0, 0.0))
quat0_l12 = [Rotations.scalar(_quat0_l12); Rotations.vector(_quat0_l12)...]
q0_l12 = [p0_l12; quat0_l12]

p1_l12 = copy(p0_l12)
quat1_l12 = copy(quat0_l12)
q1_l12 = [p1_l12; quat1_l12]

rot_off_l11_l12 = L_multiply(conjugate(quat1_l11)) * quat1_l12

idx_contact_links = [3, 6, 9]#, 12]
@assert length(idx_contact_links) == nc

function residual(z, θ, κ)

	u1 = [0.0; 0.0; 0.0]
	u2 = [0.0; 0.0; 0.0]
	u3 = [0.0; 0.0; 0.0]
	u4 = [0.0; 0.0; 0.0]
	u5 = [0.0; 0.0; 0.0]
	u6 = [0.0; 0.0; 0.0]
	u7 = [0.0; 0.0; 0.0]
	u8 = [0.0; 0.0; 0.0]
	u9 = [0.0; 0.0; 0.0]
	u10 = [0.0; 0.0; 0.0]
	u11 = [0.0; 0.0; 0.0]
	u12 = [0.0; 0.0; 0.0]

	# floating base
	q2_fb = z[1:nq]

	# link s
	q2_links = [z[nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	# joint forces
	f1 = [z[nq + N * nq + (i - 1) * 3 .+ (1:3)] for i = 1:nf]

	# joint moments
	w1 = [z[nq + N * nq + nf * 3 + (i - 1) * 2 .+ (1:2)] for i = 1:nr]

	# impact impulses (and associated slacks)
	γ1 = z[nq + N * nq + nf * 3 + nr * 2 .+ (1:nc)]
	s1 = z[nq + N * nq + nf * 3 + nr * 2 + nc .+ (1:nc)]
	b1 = z[nq + N * nq + nf * 3 + nr * 2 + 2 * nc .+ (1:(nc * nb))]
	η1 = z[nq + N * nq + nf * 3 + nr * 2 + 2 * nc + nc * nb .+ (1:(nc * nb))]
	ψ1 = z[nq + N * nq + nf * 3 + nr * 2 + 2 * nc + 2 * nc * nb .+ (1:nc)]
	s2 = z[nq + N * nq + nf * 3 + nr * 2 + 2 * nc + 2 * nc * nb + nc .+ (1:nc)]

	# f_pin = z[nq + N * nq + nf * 3 + nr * 2 + 2 * nc + 2 * nc * nb + 2 * nc .+ (1:3)]
	# w_pin = z[nq + N * nq + nf * 3 + nr * 2 + 2 * nc + 2 * nc * nb + 2 * nc + 3 .+ (1:3)]

	# previous configurations
	q0_fb = θ[1:nq]
	q0_links = [θ[nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	q1_fb = θ[nq + N * nq .+ (1:nq)]
	q1_links = [θ[nq + N * nq + nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	# time step
	h = θ[nq + N * nq + nq + N * nq .+ (1:1)]

	## contact
	#
	# signed distance
	ϕ1 = [k_func(link.kinematics[2], q2_links[i])[3] for i in idx_contact_links]

	# contact forces
	λ1 = [[friction_map * b1[(i - 1) * nb .+ (1:nb)]; γ1[i]] for i = 1:nc]
	# λ1 = [[zeros(2); γ1[i]] for i = 1:nc]

	# tangential velocity
	v = [(q2_links[i][1:3] - q1_links[i][1:3]) / h for i in idx_contact_links]
	ω = [ω_finite_difference(q1_links[i][4:7], q2_links[i][4:7], h) for i in idx_contact_links]

	P = [∇k_func(link.kinematics[2], q2_links[i])[1:2, :] for i in idx_contact_links]
	vT = [P[i] * [v[i]; ω[i]] for i = 1:nc]
	vT_stack = vcat([[vT[i]; -vT[i]] for i = 1:nc]...)

	ψ_stack = vcat([ψ1[i] * ones(nb) for i = 1:nc]...)

	# angular velocities | joint friction
	μ_joint = 1.0

	ω1 = rotation_angular_velocity(get_quaternion(q1_fb), get_quaternion(q1_links[1]), get_quaternion(q2_fb), get_quaternion(q2_links[1]), rot_off_fb_l1, h)
	ω2 = rotation_angular_velocity(get_quaternion(q1_links[1]), get_quaternion(q1_links[2]), get_quaternion(q2_links[1]), get_quaternion(q2_links[2]), rot_off_l1_l2, h)
	ω3 = rotation_angular_velocity(get_quaternion(q1_links[2]), get_quaternion(q1_links[3]), get_quaternion(q2_links[2]), get_quaternion(q2_links[3]), rot_off_l2_l3, h)

	ω4 = rotation_angular_velocity(get_quaternion(q1_fb), get_quaternion(q1_links[4]), get_quaternion(q2_fb), get_quaternion(q2_links[4]), rot_off_fb_l4, h)
	ω5 = rotation_angular_velocity(get_quaternion(q1_links[4]), get_quaternion(q1_links[5]), get_quaternion(q2_links[4]), get_quaternion(q2_links[5]), rot_off_l4_l5, h)
	ω6 = rotation_angular_velocity(get_quaternion(q1_links[5]), get_quaternion(q1_links[6]), get_quaternion(q2_links[5]), get_quaternion(q2_links[6]), rot_off_l5_l6, h)

	ω7 = rotation_angular_velocity(get_quaternion(q1_fb), get_quaternion(q1_links[7]), get_quaternion(q2_fb), get_quaternion(q2_links[7]), rot_off_fb_l7, h)
	ω8 = rotation_angular_velocity(get_quaternion(q1_links[7]), get_quaternion(q1_links[8]), get_quaternion(q2_links[7]), get_quaternion(q2_links[8]), rot_off_l7_l8, h)
	ω9 = rotation_angular_velocity(get_quaternion(q1_links[8]), get_quaternion(q1_links[9]), get_quaternion(q2_links[8]), get_quaternion(q2_links[9]), rot_off_l8_l9, h)

	# ω10 = rotation_angular_velocity(get_quaternion(q1_fb), get_quaternion(q1_links[10]), get_quaternion(q2_fb), get_quaternion(q2_links[10]), rot_off_fb_l10, h)
	# ω11 = rotation_angular_velocity(get_quaternion(q1_links[10]), get_quaternion(q1_links[11]), get_quaternion(q2_links[10]), get_quaternion(q2_links[11]), rot_off_l10_l11, h)
	# ω12 = rotation_angular_velocity(get_quaternion(q1_links[11]), get_quaternion(q1_links[12]), get_quaternion(q2_links[11]), get_quaternion(q2_links[12]), rot_off_l11_l12, h)

	μ_joint = 1.0

	joint_1_friction = -μ_joint * h[1] * ω1
	joint_2_friction = -μ_joint * h[1] * ω2
	joint_3_friction = -μ_joint * h[1] * ω3

	joint_4_friction = -μ_joint * h[1] * ω4
	joint_5_friction = -μ_joint * h[1] * ω5
	joint_6_friction = -μ_joint * h[1] * ω6

	joint_7_friction = -μ_joint * h[1] * ω7
	joint_8_friction = -μ_joint * h[1] * ω8
	joint_9_friction = -μ_joint * h[1] * ω9

	# joint_10_friction = -μ_joint * h[1] * ω10
	# joint_11_friction = -μ_joint * h[1] * ω11
	# joint_12_friction = -μ_joint * h[1] * ω12


	[
	 # floating base
	 d_fb_func(h, q0_fb, q1_fb, q2_fb,
	 	-transpose(∇k_func(floating_base.kinematics[1], q2_fb)) * f1[1]
			-transpose(∇k_func(floating_base.kinematics[2], q2_fb)) * f1[4]
			-transpose(∇k_func(floating_base.kinematics[3], q2_fb)) * f1[7],
			# -transpose(∇k_func(floating_base.kinematics[4], q2_fb)) * f1[10],
		+ transpose(dra_parent_func(get_quaternion(q2_fb), get_quaternion(q2_links[1]), rot_off_fb_l1)) * (x_axis_mask * w1[1] + x_axis_input * joint_1_friction)
   			+ transpose(dra_parent_func(get_quaternion(q2_fb), get_quaternion(q2_links[4]), rot_off_fb_l4)) * (x_axis_mask * w1[4] + x_axis_input * joint_4_friction)
			+ transpose(dra_parent_func(get_quaternion(q2_fb), get_quaternion(q2_links[7]), rot_off_fb_l7)) * (x_axis_mask * w1[7] + x_axis_input * joint_7_friction)
			# + transpose(dra_parent_func(get_quaternion(q2_fb), get_quaternion(q2_links[10]), rot_off_fb_l10)) * (x_axis_mask * w1[10] + x_axis_input * joint_10_friction)
			-u1 - u4 - u7 - u10); # floating base dynamics

	 # leg 1
	 d_link_func(h, q0_links[1], q1_links[1], q2_links[1],
	 	transpose(∇k_func(link.kinematics[2], q2_links[1])) * f1[1]
			- transpose(∇k_func(link.kinematics[1], q2_links[1])) * f1[2],
		transpose(dra_child_func(get_quaternion(q2_links[1]), get_quaternion(q2_fb), rot_off_fb_l1)) * (x_axis_mask * w1[1] + x_axis_input * joint_1_friction)
			+ transpose(dra_parent_func(get_quaternion(q2_links[1]), get_quaternion(q2_links[2]), rot_off_l1_l2)) * (z_axis_mask * w1[2] + z_axis_input * joint_2_friction)
			+ u1 - u2); # link 1 dynamics
	 d_link_func(h, q0_links[2], q1_links[2], q2_links[2],
	 	transpose(∇k_func(link.kinematics[1], q2_links[2])) * f1[2]
			- transpose(∇k_func(link.kinematics[2], q2_links[2])) * f1[3],
	 	transpose(dra_child_func(get_quaternion(q2_links[2]), get_quaternion(q2_links[1]), rot_off_l1_l2)) * (z_axis_mask * w1[2] + z_axis_input * joint_2_friction)
			+ transpose(dra_parent_func(get_quaternion(q2_links[2]), get_quaternion(q2_links[3]), rot_off_l2_l3)) * (y_axis_mask * w1[3] + y_axis_input * joint_3_friction)
			+ u2 - u3); # link 2 dynamics
	 d_link_func(h, q0_links[3], q1_links[3], q2_links[3],
	 	transpose(∇k_func(link.kinematics[1], q2_links[3])) * f1[3]
			+ transpose(∇k_func(link.kinematics[2], q2_links[3])) * λ1[1],
		transpose(dra_child_func(get_quaternion(q2_links[3]), get_quaternion(q2_links[2]), rot_off_l2_l3)) * (y_axis_mask * w1[3] + y_axis_input * joint_3_friction)
			+ u3); # link 3 dynamics

	# leg 2
	d_link_func(h, q0_links[4], q1_links[4], q2_links[4],
   	 	transpose(∇k_func(link.kinematics[2], q2_links[4])) * f1[4]
			- transpose(∇k_func(link.kinematics[1], q2_links[4])) * f1[5],
   		transpose(dra_child_func(get_quaternion(q2_links[4]), get_quaternion(q2_fb), rot_off_fb_l4)) * (x_axis_mask * w1[4] + x_axis_input * joint_4_friction)
   			+ transpose(dra_parent_func(get_quaternion(q2_links[4]), get_quaternion(q2_links[5]), rot_off_l4_l5)) * (z_axis_mask * w1[5] + z_axis_input * joint_5_friction)
   			+ u4 - u5); # link 4 dynamics
   	 d_link_func(h, q0_links[5], q1_links[5], q2_links[5],
   	 	transpose(∇k_func(link.kinematics[1], q2_links[5])) * f1[5]
			- transpose(∇k_func(link.kinematics[2], q2_links[5])) * f1[6],
   	 	transpose(dra_child_func(get_quaternion(q2_links[5]), get_quaternion(q2_links[4]), rot_off_l4_l5)) * (z_axis_mask * w1[5] + z_axis_input * joint_5_friction)
   			+ transpose(dra_parent_func(get_quaternion(q2_links[5]), get_quaternion(q2_links[6]), rot_off_l5_l6)) * (y_axis_mask * w1[6] + y_axis_input * joint_6_friction)
   			+ u5 - u6); # link 5 dynamics
   	 d_link_func(h, q0_links[6], q1_links[6], q2_links[6],
   	 	transpose(∇k_func(link.kinematics[1], q2_links[6])) * f1[6]
   			+ transpose(∇k_func(link.kinematics[2], q2_links[6])) * λ1[2],
   		transpose(dra_child_func(get_quaternion(q2_links[6]), get_quaternion(q2_links[5]), rot_off_l5_l6)) * (y_axis_mask * w1[6] + y_axis_input * joint_6_friction)
   			+ u6); # link 6 dynamics

	 # leg 3
	 d_link_func(h, q0_links[7], q1_links[7], q2_links[7],
	 	transpose(∇k_func(link.kinematics[2], q2_links[7])) * f1[7]
			- transpose(∇k_func(link.kinematics[1], q2_links[7])) * f1[8],
		transpose(dra_child_func(get_quaternion(q2_links[7]), get_quaternion(q2_fb), rot_off_fb_l7)) * (x_axis_mask * w1[7] + x_axis_input * joint_7_friction)
			+ transpose(dra_parent_func(get_quaternion(q2_links[7]), get_quaternion(q2_links[8]), rot_off_l7_l8)) * (z_axis_mask * w1[8] + z_axis_input * joint_8_friction)
			+ u7 - u8); # link 7 dynamics
	 d_link_func(h, q0_links[8], q1_links[8], q2_links[8],
	 	transpose(∇k_func(link.kinematics[1], q2_links[8])) * f1[8]
			- transpose(∇k_func(link.kinematics[2], q2_links[8])) * f1[9],
	 	transpose(dra_child_func(get_quaternion(q2_links[8]), get_quaternion(q2_links[7]), rot_off_l7_l8)) * (z_axis_mask * w1[8] + z_axis_input * joint_8_friction)
			+ transpose(dra_parent_func(get_quaternion(q2_links[8]), get_quaternion(q2_links[9]), rot_off_l8_l9)) * (y_axis_mask * w1[9] + y_axis_input * joint_9_friction)
			+ u8 - u9); # link 8 dynamics
	 d_link_func(h, q0_links[9], q1_links[9], q2_links[9],
	 	transpose(∇k_func(link.kinematics[1], q2_links[9])) * f1[9]
			+ transpose(∇k_func(link.kinematics[2], q2_links[9])) * λ1[3],
		transpose(dra_child_func(get_quaternion(q2_links[9]), get_quaternion(q2_links[8]), rot_off_l8_l9)) * (y_axis_mask * w1[9] + y_axis_input * joint_9_friction)
			+ u9); # link 9 dynamics
	#
	# # leg 4
	# d_link_func(h, q0_links[10], q1_links[10], q2_links[10],
	#    transpose(∇k_func(link.kinematics[2], q2_links[10])) * f1[10]
	# 	   - transpose(∇k_func(link.kinematics[1], q2_links[10])) * f1[11],
	#    transpose(dra_child_func(get_quaternion(q2_links[10]), get_quaternion(q2_fb), rot_off_fb_l10)) * (x_axis_mask * w1[10] + x_axis_input * joint_10_friction)
	# 	   + transpose(dra_parent_func(get_quaternion(q2_links[10]), get_quaternion(q2_links[11]), rot_off_l10_l11)) * (z_axis_mask * w1[11] + z_axis_input * joint_11_friction)
	# 	   + u10 - u11); # link 10 dynamics
	# d_link_func(h, q0_links[11], q1_links[11], q2_links[11],
	#    transpose(∇k_func(link.kinematics[1], q2_links[11])) * f1[11]
	# 	   - transpose(∇k_func(link.kinematics[2], q2_links[11])) * f1[12],
	#    transpose(dra_child_func(get_quaternion(q2_links[11]), get_quaternion(q2_links[10]), rot_off_l10_l11)) * (z_axis_mask * w1[11] + z_axis_input * joint_11_friction)
	# 	   + transpose(dra_parent_func(get_quaternion(q2_links[11]), get_quaternion(q2_links[12]), rot_off_l11_l12)) * (y_axis_mask * w1[12] + y_axis_input * joint_12_friction)
	# 	   + u11 - u12); # link 8 dynamics
	# d_link_func(h, q0_links[12], q1_links[12], q2_links[12],
	#    transpose(∇k_func(link.kinematics[1], q2_links[12])) * f1[12]
	# 	   + transpose(∇k_func(link.kinematics[2], q2_links[12])) * λ1[4],
	#    transpose(dra_child_func(get_quaternion(q2_links[12]), get_quaternion(q2_links[11]), rot_off_l11_l12)) * (y_axis_mask * w1[12] + y_axis_input * joint_12_friction)
	# 	   + u12); # link 12 dynamics

	 k_func(floating_base.kinematics[1], q2_fb) - k_func(link.kinematics[2], q2_links[1]); # body to link 1
	 k_func(link.kinematics[1], q2_links[1]) - k_func(link.kinematics[1], q2_links[2]); # link 1 to link 2
	 k_func(link.kinematics[2], q2_links[2]) - k_func(link.kinematics[1], q2_links[3]); # link 2 to link 3
	 k_func(floating_base.kinematics[2], q2_fb) - k_func(link.kinematics[2], q2_links[4]); # body to link 4
	 k_func(link.kinematics[1], q2_links[4]) - k_func(link.kinematics[1], q2_links[5]); # link 4 to link 5
	 k_func(link.kinematics[2], q2_links[5]) - k_func(link.kinematics[1], q2_links[6]); # link 5 to link 6
	 k_func(floating_base.kinematics[3], q2_fb) - k_func(link.kinematics[2], q2_links[7]); # body to link 7
	 k_func(link.kinematics[1], q2_links[7]) - k_func(link.kinematics[1], q2_links[8]); # link 7 to link 8
	 k_func(link.kinematics[2], q2_links[8]) - k_func(link.kinematics[1], q2_links[9]); # link 8 to link 9
	 # k_func(floating_base.kinematics[4], q2_fb) - k_func(link.kinematics[2], q2_links[10]); # body to link 10
	 # k_func(link.kinematics[1], q2_links[10]) - k_func(link.kinematics[1], q2_links[11]); # link 10 to link 11
	 # k_func(link.kinematics[2], q2_links[11]) - k_func(link.kinematics[1], q2_links[12]); # link 11 to link 12

	 transpose(x_axis_mask) * ra_func(get_quaternion(q2_fb), get_quaternion(q2_links[1]), rot_off_fb_l1);
	 transpose(z_axis_mask) * ra_func(get_quaternion(q2_links[1]), get_quaternion(q2_links[2]), rot_off_l1_l2);
	 transpose(y_axis_mask) * ra_func(get_quaternion(q2_links[2]), get_quaternion(q2_links[3]), rot_off_l2_l3);
	 transpose(x_axis_mask) * ra_func(get_quaternion(q2_fb), get_quaternion(q2_links[4]), rot_off_fb_l4);
	 transpose(z_axis_mask) * ra_func(get_quaternion(q2_links[4]), get_quaternion(q2_links[5]), rot_off_l4_l5);
	 transpose(y_axis_mask) * ra_func(get_quaternion(q2_links[5]), get_quaternion(q2_links[6]), rot_off_l5_l6);
	 transpose(x_axis_mask) * ra_func(get_quaternion(q2_fb), get_quaternion(q2_links[7]), rot_off_fb_l7);
	 transpose(z_axis_mask) * ra_func(get_quaternion(q2_links[7]), get_quaternion(q2_links[8]), rot_off_l7_l8);
	 transpose(y_axis_mask) * ra_func(get_quaternion(q2_links[8]), get_quaternion(q2_links[9]), rot_off_l8_l9);
	 # transpose(x_axis_mask) * ra_func(get_quaternion(q2_fb), get_quaternion(q2_links[10]), rot_off_fb_l10);
	 # transpose(z_axis_mask) * ra_func(get_quaternion(q2_links[10]), get_quaternion(q2_links[11]), rot_off_l10_l11);
	 # transpose(y_axis_mask) * ra_func(get_quaternion(q2_links[11]), get_quaternion(q2_links[12]), rot_off_l11_l12);

	 s1 .- ϕ1;
	 η1 .- vT_stack .- ψ_stack;
	 s2 .- (μ_world .* γ1 .- E * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ

	 ]
end
nz = nq + N * nq + nf * 3 + nr * 2 + n_contact
nθ = 2 * nq + N * (2 * nq) + 1
@variables z[1:nz], θ[1:nθ], κ[1:1]

tmp = Array(Diagonal(ones(nf * 3 + nr * 2 + n_contact)))
function Gz_func(q)
	q_fb = q[1:nq]
	q_links = [q[nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	cat(G_func(q_fb), [G_func(q_links[i]) for i = 1:N]..., tmp, dims=(1, 2))
end

r = vec(residual(z, θ, κ))
∇r = Symbolics.jacobian(r, z) * Gz_func(z)

r_func = eval(Symbolics.build_function(r, z, θ, κ)[2])
∇r_func = eval(Symbolics.build_function(∇r, z, θ)[2])

rz = similar(∇r, Float64)

rq_space = rn_quaternion_space(nz - (1 + N), x -> Gz_func(x),
	vcat(vcat([(i - 1) * nq .+ (1:3) for i = 1:(1 + N)]...), collect((nq * (1 + N) .+ (1:(3 * nf + 2 * nr + n_contact))))),
	vcat(vcat([(i - 1) * (nq - 1) .+ (1:3) for i = 1:(1 + N)]...), collect(((nq - 1) * (1 + N) .+ (1:(3 * nf + 2 * nr + n_contact))))),
	[collect((i - 1) * nq + 3 .+ (1:4)) for i = 1:(1 + N)],
	[collect((i - 1) * (nq - 1) + 3 .+ (1:3)) for i = 1:(1 + N)])


# options
opts = ContactControl.InteriorPointOptions(diff_sol = false, r_tol = 1.0e-8, κ_tol = 1.0e-6)

r_func(zeros(nz - (1 + N)), ones(nz), ones(nθ), 1.0)
∇r_func(rz, ones(nz), ones(nθ))

h = 0.01
θ0 = [q0_fb;
	  q0_l1; q0_l2; q0_l3;
	  q0_l4; q0_l5; q0_l6;
	  q0_l7; q0_l8; q0_l9;
	  # q0_l10; q0_l11; q0_l12;
	  q1_fb;
	  q1_l1; q1_l2; q1_l3;
	  q1_l4; q1_l5; q1_l6;
	  q1_l7; q1_l8; q1_l9;
	  # q1_l10; q1_l11; q1_l12;
	  h]

z0 = copy([q1_fb;
		   q1_l1; q1_l2; q1_l3;
		   q1_l4; q1_l5; q1_l6;
		   q1_l7; q1_l8; q1_l9;
		   # q1_l10; q1_l11; q1_l12;
		   zeros(nf * 3 + nr * 2); 0.1 * ones(n_contact)])

# _r = zeros(36)
# r_func(_r, z0, θ0, 0.0)
# _r
# solver
idx_ineq = collect(nq + N * nq + nf * 3 + nr * 2 .+ (1:(n_contact)))
ip = ContactControl.interior_point(z0, θ0,
	s = rq_space,
	idx_ineq = idx_ineq,
	r! = r_func, rz! = ∇r_func,
	rz = rz,
	opts = opts)

# solve
T = 50
q_hist = [[q0_fb;
		   q0_l1; q0_l2; q0_l3;
		   q0_l4; q0_l5; q0_l6;
		   q0_l7; q0_l8; q0_l9;
		   # q0_l10; q0_l11; q0_l12
		   ],
		  [q1_fb;
		   q1_l1; q1_l2; q1_l3;
		   q1_l4; q1_l5; q1_l6;
		   q1_l7; q1_l8; q1_l9;
		   # q1_l10; q1_l11; q1_l12
		   ]]

ip.opts.κ_init = 0.1
for t = 1:T-2
	println("t = $t")
	ip.θ .= [q_hist[end-1]; q_hist[end]; h]
	ip.z .= copy([q_hist[end]; 0.001 * randn(nf * 3 + nr * 2); 0.1 * ones(n_contact)])
	status = ContactControl.interior_point_solve!(ip)
	if !status
		println("t: ", t)
		println("res norm: ", norm(ip.r, Inf))
	end
	push!(q_hist, ip.z[1:((1 + N) * nq)])
end

function visualize!(vis, fb::Body, links::Vector{Body}, q;
        Δt = 0.1, r_link = 0.025)

	default_background!(vis)

	N = length(links)

    setobject!(vis["fb"], GeometryBasics.Rect(Vec(-1.0 * 0.5 * fb_length,
		-1.0 * 0.5 * fb_width,
		-1.0 * 0.5 * fb_height),
		Vec(2.0 * 0.5 * fb_length, 2.0 * 0.5 * fb_width, 2.0 * 0.5 * fb_height)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 0.5)))

	for i = 1:N
		setobject!(vis["l$i"], GeometryBasics.Rect(Vec(-r_link,
			-r_link,
			-1.0 * 0.5 * link_length),
			Vec(2.0 * r_link, 2.0 * r_link, 2.0 * 0.5 * link_length)),
			MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 0.25)))
	end

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["fb"],
				compose(Translation(q[t][1:3]...), LinearMap(UnitQuaternion(q[t][4:7]...))))
			for i = 1:N
				settransform!(vis["l$i"],
					compose(Translation(q[t][nq + (i - 1) * nq .+ (1:3)]...), LinearMap(UnitQuaternion(q[t][nq + (i - 1) * nq .+ (4:7)]...))))
			end
        end
    end
    MeshCat.setanimation!(vis, anim)
end

vis = Visualizer()
render(vis)
visualize!(vis, floating_base, links, q_hist, Δt = h)
