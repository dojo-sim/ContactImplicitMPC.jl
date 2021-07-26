struct Body
	mass
	inertia
	length
	gravity
end

function rotation_axes(q1, q2)
	return multiply(conjugate(q1), q2)[2:4]
end

function get_quaternion(q)
	return q[4:7]
end

function G_func(::Body, q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

body = Body(Diagonal(1.0 * ones(3)),
			Diagonal([0.1, 0.1, 0.1]),
			0.5,
			[0.0; 0.0; 0.0 * 9.81])
nq = 7

@variables q[1:nq]

function kinematics1(body, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	# body end effector
	pe1 = [0.0; 0.0; body.length]

	k1 = p + quaternion_rotation_matrix(quat) * pe1

	return k1
end

function kinematics2(body, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	# body end effector
	pe2 = [0.0; 0.0; -body.length]

	k2 = p + quaternion_rotation_matrix(quat) * pe2

	return k2
end

k1 = kinematics1(body, q)
k2 = kinematics2(body, q)

k1_func = eval(Symbolics.build_function(k1, q)[1])
k2_func = eval(Symbolics.build_function(k2, q)[1])

∇k1 = Symbolics.jacobian(k1, q) * G_func(body, q)
∇k2 = Symbolics.jacobian(k2, q) * G_func(body, q)

∇k1_func = eval(Symbolics.build_function(∇k1, q)[1])
∇k2_func = eval(Symbolics.build_function(∇k2, q)[1])


function dynamics(model::Body, h, q0, q1, q2, f1, f2, τ1)

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
	d_angular = (model.inertia * ω2 * sqrt(4.0 / h[1]^2.0 - transpose(ω2) * ω2)
		+ cross(ω2, model.inertia * ω2)
		- model.inertia * ω1 * sqrt(4.0 / h[1]^2.0 - transpose(ω1) * ω1)
		+ cross(ω1, model.inertia * ω1) - 2.0 * transpose(quaternion_rotation_matrix(quat2)) * τ1)

	return [d_linear; d_angular] + transpose(∇k1_func(q2)) * f1 + transpose(∇k2_func(q2)) * f2
end

function G_func(::Body, q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

@variables q0[1:nq], q1[1:nq], q2[1:nq], f1[1:3], f2[1:3], τ1[1:3], h[1:1]

d = dynamics(body, h, q0, q1, q2, f1, f2, τ1)
G = G_func(body, q2)
∇d = Symbolics.jacobian(d, q2) * G

d_func = eval(Symbolics.build_function(d, h, q0, q1, q2, f1, f2, τ1)[1])
∇d_func = eval(Symbolics.build_function(∇d, h, q0, q1, q2, f1, f2, τ1)[1])

d_func([1.0], ones(7), ones(7), ones(7), ones(3), ones(3), ones(3))

x_mask_axis = [0.0 0.0;
               1.0 0.0;
			   0.0 1.0]

function residual(z, θ, κ)
	q2_link1 = z[1:7]
	q2_link2 = z[7 .+ (1:7)]

	f_world_to_link1 = z[2 * 7 .+ (1:3)]
	f_link1_to_link2 = z[2 * 7 + 3 .+ (1:3)]

	τ_world_to_link1 = z[2 * 7 + 2 * 3 .+ (1:2)]
	# τ_link1_to_link2 = z[2 * 7 + 2 * 3 + 2 .+ (1:2)]

	q0_link1 = θ[1:7]
	q0_link2 = θ[7 .+ (1:7)]

	q1_link1 = θ[2 * 7 .+ (1:7)]
	q1_link2 = θ[3 * 7 .+ (1:7)]

	f_link2 = θ[4 * 7 .+ (1:3)]
	τ_link1 = θ[4 * 7 + 3 .+ (1:3)]
	τ_link2 = θ[4 * 7 +  2 * 3 .+ (1:3)]

	rot_world_l1 = rotation_axes(get_quaternion(q2_link1), [1.0; 0.0; 0.0; 0.0])
	rot_world_l1 = rotation_axes(get_quaternion(q2_link1), [1.0; 0.0; 0.0; 0.0])

	h = θ[4 * 7 + 3 * 3 .+ (1:1)]

	[
	 d_func(h, q0_link1, q1_link1, q2_link1, f_world_to_link1, f_link1_to_link2, τ_link1); # link 1 dynamics
	 d_func(h, q0_link2, q1_link2, q2_link2, -f_link1_to_link2, f_link2, τ_link2) # link 2 dynamics
	 k2_func(q2_link1); # link 1 to world
	 k1_func(q2_link1) - k2_func(q2_link2)
	 ]

	# [d_func(h, q0_link1, q1_link1, q2_link1, f_link1_to_link2, f_world_to_link1, τ_link1); # link 1 dynamics
	#  d_func(h, q0_link2, q1_link2, q2_link2, f_link2, -f_link1_to_link2, τ_link2) # link 2 dynamics
	#  k2_func(q2_link1); # link 1 to world
	#  k1_func(q2_link1) - k2_func(q2_link2) # link 1 to link 2
	#  ]
end

num_bodies = 2
nz = num_bodies * nq + 3 * 2
nθ = num_bodies * (2 * nq) + 3 * 3 + 1
@variables z[1:nz], θ[1:nθ], κ[1:1]

tmp = Array(Diagonal(ones(6)))
function Gz_func(body::Body, q)
	q_link1 = q[1:nq]
	q_link2 = q[nq .+ (1:nq)]

	cat(G_func(body, q_link1), G_func(body, q_link2),  tmp, dims=(1, 2))
end

r = residual(z, θ, κ)
∇r = Symbolics.jacobian(r, z) * Gz_func(body, z)

r_func = eval(Symbolics.build_function(r, z, θ, κ)[2])
∇r_func = eval(Symbolics.build_function(∇r, z, θ)[2])

rz = similar(∇r, Float64)

# options
opts = ContactControl.InteriorPointOptions(diff_sol = false)

r_func(zeros(nz-2), ones(nz), ones(nθ), 1.0)
∇r_func(rz, ones(nz), ones(nθ))

# Rn + quaternion space
rq_space = rn_quaternion_space(nz - 2, x -> Gz_func(body, x),
			collect([(1:3)..., (8:10)..., (15:20)...]),
			collect([(1:3)..., (7:9)..., (13:18)...]),
			[collect(4:7), collect(11:14)],
			[collect(4:6), collect(10:12)])

p0_link1 = [0.0; body.length; 0.0]
# p0_link1 = [0.0; 0.0; body.length]
_quat0_link1 = UnitQuaternion(AngleAxis(-0.5 * π, 1.0, 0.0, 0.0))
# _quat0_link1 = one(UnitQuaternion)
quat0_link1 = [Rotations.scalar(_quat0_link1); Rotations.vector(_quat0_link1)...]
q0_link1 = [p0_link1; quat0_link1]

p1_link1 = [0.0; body.length; 0.0]
# p0_link1 = [0.0; 0.0; body.length]
quat1_link1 = copy(quat0_link1)
q1_link1 = [p1_link1; quat1_link1]

p0_link2 = [0.0; 3 * body.length; 0.0]
# p0_link2 = [0.0; 0.0; 3 * body.length]
_quat0_link2 = UnitQuaternion(AngleAxis(-0.5 * π, 1.0, 0.0, 0.0))
# _quat0_link2 = one(UnitQuaternion)
quat0_link2 = [Rotations.scalar(_quat0_link2); Rotations.vector(_quat0_link2)...]
q0_link2 = [p0_link2; quat0_link2]

p1_link2 = [0.0; 3 * body.length; 0.0]
# p1_link2 = [0.0; 0.0; 3 * body.length]
quat1_link2 = copy(quat0_link2)
q1_link2 = [p1_link2; quat1_link2]

h = 0.1
f_link2 = zeros(3)
τ_link1 = [0.0; 0.0; 0.0]
τ_link2 = [0.0; 0.0; 0.0]
θ0 = [q0_link1; q0_link2; q1_link1; q1_link2; f_link2; τ_link1; τ_link2; h]
z0 = copy([q1_link1; q1_link2; 0.0 * randn(6)])

# solver
ip = ContactControl.interior_point(z0, θ0,
	s = rq_space,
	idx_ineq = collect(1:0),
	r! = r_func, rz! = ∇r_func,
	rz = rz,
	opts = opts)

# solve
T = 20
q_hist = [[q0_link1; q0_link2], [q1_link1; q1_link2]]

for t = 1:T-2
	ip.θ .= [q_hist[end-1]; q_hist[end]; f_link2; τ_link1; τ_link2; h]
	ip.z .= copy([q_hist[end]; 0.0 * randn(6)])
	status = ContactControl.interior_point_solve!(ip)
	push!(q_hist, ip.z[1:(2 * nq)])
end

function visualize!(vis, model::Body, q;
        Δt = 0.1)

	default_background!(vis)

	r = model.length

    setobject!(vis["link1"], GeometryBasics.Rect(Vec(-0.1 * r,
		-0.1 * r,
		-1.0 * r),
		Vec(0.2 * r, 0.2 * r, 2.0 * r)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	setobject!(vis["link2"], GeometryBasics.Rect(Vec(-0.1 * r,
		-0.1 * r,
		-1.0 * r),
		Vec(0.2 * r, 0.2 * r, 2.0 * r)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["link1"],
				compose(Translation(q[t][1:3]...), LinearMap(UnitQuaternion(q[t][4:7]...))))
			settransform!(vis["link2"],
				compose(Translation(q[t][7 .+ (1:3)]...), LinearMap(UnitQuaternion(q[t][7 .+ (4:7)]...))))
        end
    end
    MeshCat.setanimation!(vis, anim)
end

vis = Visualizer()
render(vis)
visualize!(vis, body, q_hist, Δt = h)
