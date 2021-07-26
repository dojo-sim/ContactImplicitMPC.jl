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

floating_base = Body(Diagonal(1.0 * ones(3)),
			Diagonal([0.1, 0.1, 0.1]),
			[0.5, 0.25, 0.1],
			[0.0; 0.0; 0.0 * 9.81])

link = Body(Diagonal(1.0 * ones(3)),
			Diagonal([0.1, 0.1, 0.1]),
			0.25,
			[0.0; 0.0; 0.0 * 9.81])

nq = 7

@variables q[1:nq]

function kinematics1(body, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	# body end effector
	pe1 = [0.0; 0.0; 0.5 * body.length]

	k1 = p + quaternion_rotation_matrix(quat) * pe1

	return k1
end

function kinematics2(body, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	# body end effector
	pe2 = [0.0; 0.0; -0.5 * body.length]

	k2 = p + quaternion_rotation_matrix(quat) * pe2

	return k2
end

function kinematics_shoulder1(body, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	# body end effector
	pe1 = [0.5 * body.length[1]; 0.5 * body.length[2]; 0.0]

	k1 = p + quaternion_rotation_matrix(quat) * pe1

	return k1
end

k1 = kinematics1(link, q)
k2 = kinematics2(link, q)
ks1 = kinematics_shoulder1(floating_base, q)

k1_func = eval(Symbolics.build_function(k1, q)[1])
k2_func = eval(Symbolics.build_function(k2, q)[1])
ks1_func = eval(Symbolics.build_function(ks1, q)[1])

∇k1 = Symbolics.jacobian(k1, q) * G_func(link, q)
∇k2 = Symbolics.jacobian(k2, q) * G_func(link, q)
∇ks1 = Symbolics.jacobian(ks1, q) * G_func(floating_base, q)


∇k1_func = eval(Symbolics.build_function(∇k1, q)[1])
∇k2_func = eval(Symbolics.build_function(∇k2, q)[1])
∇ks1_func = eval(Symbolics.build_function(∇ks1, q)[1])


function dynamics_link(model::Body, h, q0, q1, q2, f1, f2, τ1)

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

function dynamics_base(model::Body, h, q0, q1, q2, f1, f2, f3, f4, τ1)

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

	return [d_linear; d_angular] + transpose(∇ks1_func(q2)) * f1
end

function G_func(::Body, q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

@variables q0[1:nq], q1[1:nq], q2[1:nq], f1[1:3], f2[1:3], f3[1:3], f4[1:3], τ1[1:3], h[1:1]

d_link = dynamics_link(link, h, q0, q1, q2, f1, f2, τ1)
G = G_func(link, q2)
∇d_link = Symbolics.jacobian(d_link, q2) * G

d_link_func = eval(Symbolics.build_function(d_link, h, q0, q1, q2, f1, f2, τ1)[1])
∇d_link_func = eval(Symbolics.build_function(∇d_link, h, q0, q1, q2, f1, f2, τ1)[1])

d_link_func([1.0], ones(7), ones(7), ones(7), ones(3), ones(3), ones(3))

d_fb = dynamics_base(floating_base, h, q0, q1, q2, f1, f2, f3, f4, τ1)
G = G_func(floating_base, q2)
∇d_fb = Symbolics.jacobian(d_fb, q2) * G

d_fb_func = eval(Symbolics.build_function(d_fb, h, q0, q1, q2, f1, f2, f3, f4, τ1)[1])
∇d_fb_func = eval(Symbolics.build_function(∇d_fb, h, q0, q1, q2, f1, f2, f3, f4, τ1)[1])

d_fb_func([1.0], ones(7), ones(7), ones(7), ones(3), ones(3), ones(3), ones(3), ones(3))

links = [link, link, link]
N = length(links)
nf = 3
nr = 0

x_axis_mask = [0.0 0.0;
			   1.0 0.0;
			   0.0 1.0]

y_axis_mask = [1.0 0.0;
			   0.0 0.0;
			   0.0 1.0]

# floating base
p0_fb = [0.0; 0.0; 0.0]
_quat0_fb = one(UnitQuaternion)
quat0_fb = [Rotations.scalar(_quat0_fb); Rotations.vector(_quat0_fb)...]
q0_fb = [p0_fb; quat0_fb]

p1_fb = [0.0; 0.0; 0.0]
quat1_fb = copy(quat0_fb)
q1_fb = [p1_fb; quat1_fb]

# link 1
p0_l1 = p0_fb + [0.5 * floating_base.length[1]; 0.5 * floating_base.length[2] + 0.5 * link.length; 0.0]
_quat0_l1 = UnitQuaternion(AngleAxis(-0.5 * π, 1.0, 0.0, 0.0))
quat0_l1 = [Rotations.scalar(_quat0_l1); Rotations.vector(_quat0_l1)...]
q0_l1 = [p0_l1; quat0_l1]

p1_l1 = p1_fb + [0.5 * floating_base.length[1]; 0.5 * floating_base.length[2] + 0.5 * link.length; 0.0]
quat1_l1 = copy(quat0_l1)
q1_l1 = [p1_l1; quat1_l1]

rot_off_fb_l1 = rotation_axes(get_quaternion(q1_fb), get_quaternion(q1_l1))

# link 2
p0_l2 = p0_fb + [0.5 * floating_base.length[1]; 0.5 * floating_base.length[2] + 1.0 * link.length; -0.5 * link.length]
_quat0_l2 = one(UnitQuaternion)
quat0_l2 = [Rotations.scalar(_quat0_l2); Rotations.vector(_quat0_l2)...]
q0_l2 = [p0_l2; quat0_l2]

p1_l2 = p1_fb + [0.5 * floating_base.length[1]; 0.5 * floating_base.length[2] + 1.0 * link.length; -0.5 * link.length]
quat1_l2 = copy(quat0_l2)
q1_l2 = [p1_l2; quat1_l2]

rot_off_l1_l2 = rotation_axes(get_quaternion(q1_l1), get_quaternion(q1_l2))

# link 3
p0_l3 = p0_fb + [0.5 * floating_base.length[1]; 0.5 * floating_base.length[2] + 1.0 * link.length; -1.5 * link.length ]
_quat0_l3 = one(UnitQuaternion)
quat0_l3 = [Rotations.scalar(_quat0_l3); Rotations.vector(_quat0_l3)...]
q0_l3 = [p0_l3; quat0_l3]

p1_l3 = p1_fb + [0.5 * floating_base.length[1]; 0.5 * floating_base.length[2] + 1.0 * link.length; -1.5 * link.length]
quat1_l3 = copy(quat0_l3)
q1_l3 = [p1_l3; quat1_l3]

rot_off_l2_l3 = rotation_axes(get_quaternion(q1_l2), get_quaternion(q1_l3))

function residual(z, θ, κ)
	u1 = [0.01; 0.0; 0.0]
	u2 = [0.0; 0.1; 0.0]
	u3 = [0.0; 0.0; 0.05]

	# floating base
	q2_fb = z[1:nq]

	# link 1
	q2_links = [z[nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	f1 = [z[nq + N * nq + (i - 1) * 3 .+ (1:3)] for i = 1:nf]
	w1 = [z[nq + N * nq + nf * 3 .+ (1:2)] for i = 1:nr]

	q0_fb = θ[1:nq]
	q0_links = [θ[nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	q1_fb = θ[nq + N * nq .+ (1:nq)]
	q1_links = [θ[nq + N * nq + nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	h = θ[nq + N * nq + nq + N * nq .+ (1:1)]

	rot_axes_fb_l1 = rotation_axes(get_quaternion(q2_fb), get_quaternion(q2_links[1]))
	rot_axes_l1_l2 = rotation_axes(get_quaternion(q2_links[1]), get_quaternion(q2_links[2]))
	rot_axes_l2_l3 = rotation_axes(get_quaternion(q2_links[2]), get_quaternion(q2_links[2]))

	[
	 d_fb_func(h, q0_fb, q1_fb, q2_fb, -f1[1], zeros(3), zeros(3), zeros(3), -u1);#-x_axis_mask * w1[1] - u1); # link 1 dynamics
	 d_link_func(h, q0_links[1], q1_links[1], q2_links[1], -f1[2], f1[1], u1 - u2);#x_axis_mask * w1[1]);# - y_axis_mask * w1[2] + u1 - u2); # link 2 dynamics
	 d_link_func(h, q0_links[2], q1_links[2], q2_links[2], f1[2], -f1[3], u2 - u3);#y_axis_mask * w1[2] - y_axis_mask * w1[3] + u2 - u3); # link 2 dynamics
	 d_link_func(h, q0_links[3], q1_links[3], q2_links[3], f1[3], zeros(3), u3);#y_axis_mask * w1[3] + u3); # link 2 dynamics
	 ks1_func(q2_fb) - k2_func(q2_links[1]); # body to link 1
	 k1_func(q2_links[1]) - k1_func(q2_links[2]); # link 1 to link 2
	 k2_func(q2_links[2]) - k1_func(q2_links[3]); # link 2 to link 3
	 # transpose(x_axis_mask) * (rot_axes_fb_l1 - rot_off_fb_l1);
	 # transpose(y_axis_mask) * (rot_axes_l1_l2 - rot_off_l1_l2);
	 # transpose(y_axis_mask) * (rot_axes_l2_l3 - rot_off_l2_l3);
	 ]
end

nz = nq + N * nq + nf * 3 + nr * 2
nθ = 2 * nq + N * (2 * nq) + 1
@variables z[1:nz], θ[1:nθ], κ[1:1]

tmp = Array(Diagonal(ones(nf * 3 + nr * 2)))
function Gz_func(::Body, q)
	q_fb = q[1:nq]
	q_links = [q[nq + (i - 1) * nq .+ (1:nq)] for i = 1:N]

	cat(G_func(floating_base, q_fb), [G_func(link, q_links[i]) for i = 1:N]..., tmp, dims=(1, 2))
end

r = residual(z, θ, κ)
∇r = Symbolics.jacobian(r, z) * Gz_func(body, z)

r_func = eval(Symbolics.build_function(r, z, θ, κ)[2])
∇r_func = eval(Symbolics.build_function(∇r, z, θ)[2])

rz = similar(∇r, Float64)


# Rn + quaternion space
# rq_space = rn_quaternion_space(nz - (1 + N), x -> Gz_func(body, x),
# 			collect([(1:3)..., (8:10)..., (15:17)..., (22:24)..., (28 .+ (1:3))...]),
# 			collect([(1:3)..., (7:9)..., (13:15)..., (19:21)..., (24 .+ (1:3))...]),
# 			[collect(4:7), collect(11:14), collect(18:21), collect(25:28)],
# 			[collect(4:6), collect(10:12), collect(16:18), collect(22:24)])

rq_space = rn_quaternion_space(nz - (1 + N), x -> Gz_func(body, x),
	vcat(vcat([(i - 1) * nq .+ (1:3) for i = 1:(1 + N)]...), collect((nq * (1 + N) .+ (1:(3 * nf))))),
	vcat(vcat([(i - 1) * (nq - 1) .+ (1:3) for i = 1:(1 + N)]...), collect(((nq - 1) * (1 + N) .+ (1:(3 * nf))))),
	[collect((i - 1) * nq + 3 .+ (1:4)) for i = 1:(1 + N)],
	[collect((i - 1) * (nq - 1) + 3 .+ (1:3)) for i = 1:(1 + N)])


# options
opts = ContactControl.InteriorPointOptions(diff_sol = false)

r_func(zeros(nz-(1 + N)), ones(nz), ones(nθ), 1.0)
∇r_func(rz, ones(nz), ones(nθ))


h = 0.1

θ0 = [q0_fb; q0_l1; q0_l2; q0_l3; q1_fb; q1_l1; q1_l2; q1_l3; h]
z0 = copy([q1_fb; q1_l1; q1_l2; q1_l3; zeros(nf * 3 + nr * 2)])

# solver
ip = ContactControl.interior_point(z0, θ0,
	s = rq_space,
	idx_ineq = collect(1:0),
	r! = r_func, rz! = ∇r_func,
	rz = rz,
	opts = opts)

# solve
T = 20
q_hist = [[q0_fb; q0_l1; q0_l2; q0_l3], [q1_fb; q1_l1; q1_l2; q1_l3]]

for t = 1:T-2
	ip.θ .= [q_hist[end-1]; q_hist[end]; h]
	ip.z .= copy([q_hist[end]; zeros(nf * 3 + nr * 2)])
	# ip.z .= copy(q_hist[end])
	status = ContactControl.interior_point_solve!(ip)
	push!(q_hist, ip.z[1:((1 + N) * nq)])
end

function visualize!(vis, fb::Body, links::Vector{Body}, q;
        Δt = 0.1, r_link = 0.025)

	default_background!(vis)

	N = length(links)

    setobject!(vis["fb"], GeometryBasics.Rect(Vec(-1.0 * 0.5 * fb.length[1],
		-1.0 * 0.5 * fb.length[2],
		-1.0 * 0.5 * fb.length[3]),
		Vec(2.0 * 0.5 * fb.length[1], 2.0 * 0.5 * fb.length[2], 2.0 * 0.5 * fb.length[3])),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 0.5)))

	for i = 1:N
		setobject!(vis["l$i"], GeometryBasics.Rect(Vec(-r_link,
			-r_link,
			-1.0 * 0.5 * links[i].length),
			Vec(2.0 * r_link, 2.0 * r_link, 2.0 * 0.5 * links[i].length)),
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
