struct Body
	mass
	inertia
	length
	gravity
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

	d_linear = model.mass * (vm1[1:3] - vm2[1:3]) + h[1] * model.mass * model.gravity
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

function residual(z, θ, κ)
	q2 = z[1:7]
	f2 = z[8:10]

	q0 = θ[1:7]
	q1 = θ[7 .+ (1:7)]

	f1 = θ[14 .+ (1:3)]
	# f2 = θ[17 .+ (1:3)]
	τ1 = θ[17 .+ (1:3)]

	h = θ[21:21]

	# f1 = θ[14 .+ (1:3)]
	# f2 = θ[17 .+ (1:3)]
	# τ1 = θ[20 .+ (1:3)]
	#
	# h = θ[24:24]
	# d_func(h, q0, q1, q2, f1, zeros(3), τ1)
	[d_func(h, q0, q1, q2, f1, f2, τ1);
	 k2_func(q2)]
end

num_bodies = 1
nz = num_bodies * nq + 3
nθ = 21
@variables z[1:nz], θ[1:nθ], κ[1:1]

# tmp = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
tmp = Array(Diagonal(ones(3)))
function Gz_func(body::Body, q)
	# G_func(body, q)
	cat(G_func(body, q), tmp, dims=(1, 2))
end

r = residual(z, θ, κ)
∇r = Symbolics.jacobian(r, z) * Gz_func(body, z)

r_func = eval(Symbolics.build_function(r, z, θ, κ)[2])
∇r_func = eval(Symbolics.build_function(∇r, z, θ)[2])

rz = similar(∇r, Float64)

# options
opts = ContactControl.InteriorPointOptions(diff_sol = false)

∇r_func(rz, ones(nz), ones(nθ))

# Rn + quaternion space
rq_space = rn_quaternion_space(nz-1, x -> Gz_func(body, x),
			collect([(1:3)..., (8:10)...]),
			collect([(1:3)..., (7:9)...]),
			[collect(4:7)],
			[collect(4:6)])

p0 = [0.0; 0.0; 0.0]
_quat0 = UnitQuaternion(AngleAxis(-0.5 * π, 1.0, 0.0, 0.0))
# _quat0 = one(UnitQuaternion)
quat0 = [Rotations.scalar(_quat0); Rotations.vector(_quat0)...]
q0 = [p0; quat0]

p1 = [0.0; 0.0; 0.0]
quat1 = copy(quat0)
q1 = [p1; quat1]

f1 = [0.0; 0.0; 0.0]
f2 = [0.0; 0.0; 0.0]
τ1 = [0.1; 0.0; 0.0]

h = 0.05
θ0 = [q0; q1; f1; τ1; 0.1]
z0 = copy([q1; zeros(3)])
# z0 = copy(q1)

# solver
ip = ContactControl.interior_point(z0, θ0,
	s = rq_space,
	idx_ineq = collect(1:0),
	r! = r_func, rz! = ∇r_func,
	rz = rz,
	opts = opts)

# solve
T = 100
q_hist = [q0, q1]

for t = 1:T-2
	ip.θ .= [q_hist[end-1]; q_hist[end]; f1; τ1; h]
	ip.z .= copy([q_hist[end]; zeros(3)])
	# ip.z .= copy(q_hist[end])
	status = ContactControl.interior_point_solve!(ip)
	push!(q_hist, ip.z[1:nq])
end

function visualize!(vis, model::Body, q;
        Δt = 0.1)

	default_background!(vis)

	r = model.length

    setobject!(vis["box"], GeometryBasics.Rect(Vec(-0.1 * r,
		-0.1 * r,
		-1.0 * r),
		Vec(0.2 * r, 0.2 * r, 2.0 * r)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["box"],
				compose(Translation(q[t][1:3]...), LinearMap(UnitQuaternion(q[t][4:7]...))))
        end
    end
    MeshCat.setanimation!(vis, anim)
end

vis = Visualizer()
render(vis)
visualize!(vis, body, q_hist, Δt = h)
