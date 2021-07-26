struct Body
	mass
	inertia
	kinematics
	gravity
end

function get_quaternion(q)
	return q[4:7]
end

function rotation_axes(q1, q2)
	return multiply(conjugate(q1), q2)[2:4]
end


function G_func(q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

link = Body(Diagonal(1.0 * ones(3)),
			Diagonal([0.1, 0.1, 0.1]),
			[[0.0; 0.0; 0.5], [0.0; 0.0; -0.5]],
			[0.0; 0.0; 0.0 * 9.81])
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

@variables quat0[1:4], quat1[1:4]
ra = rotation_axes(quat0, quat1)
dra = Symbolics.jacobian(ra, quat0) * attitude_jacobian(quat0)

ra_func = eval(Symbolics.build_function(ra, quat0, quat1)[1])
dra_func = eval(Symbolics.build_function(dra, quat0, quat1)[1])

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

	d_linear = model.mass * (vm1[1:3] - vm2[1:3]) + h[1] * model.mass * model.gravity
	d_angular = (model.inertia * ω2 * sqrt(4.0 / h[1]^2.0 - transpose(ω2) * ω2)
		+ cross(ω2, model.inertia * ω2)
		- model.inertia * ω1 * sqrt(4.0 / h[1]^2.0 - transpose(ω1) * ω1)
		+ cross(ω1, model.inertia * ω1))

	d_angular .-= 2.0 * τ

	return [d_linear; d_angular] + f
end

@variables q0[1:nq], q1[1:nq], q2[1:nq], f1[1:6], τ1[1:3], h[1:1]

d = dynamics(link, h, q0, q1, q2, f1, τ1)
G = G_func(q2)
∇d = Symbolics.jacobian(d, q2) * G

d_func = eval(Symbolics.build_function(d, h, q0, q1, q2, f1, τ1)[1])
∇d_func = eval(Symbolics.build_function(∇d, h, q0, q1, q2, f1, τ1)[1])

d_func([1.0], ones(7), ones(7), ones(7), ones(6), ones(3))

x_axis_mask = [0.0 0.0;
               1.0 0.0;
			   0.0 1.0]

function residual(z, θ, κ)
	q2 = z[1:nq]
	_f = z[8:10]
	_τ = z[11:12]

	quat2 = q2[4:7]

	q0 = θ[1:7]
	q1 = θ[7 .+ (1:7)]

	h = θ[15:15]

	quat_world = [1.0; 0.0; 0.0; 0.0]
	ra = ra_func(quat2, quat_world)

	f = transpose(∇k_func(link.kinematics[2], q2)) * _f
	τ = transpose(dra_func(quat2, quat_world)) * x_axis_mask * _τ
	u = transpose(quaternion_rotation_matrix(quat2)) * [0.5; 0.0; 0.5]
	[
	 d_func(h, q0, q1, q2, f, τ + u); # dynamics
	 k_func(link.kinematics[2], q2); # position constraint
	 transpose(x_axis_mask) * ra_func(quat2, quat_world)
	 ]
end

num_bodies = 1
nz = num_bodies * nq + 3 + 2
nθ = 15
@variables z[1:nz], θ[1:nθ], κ[1:1]

tmp = Array(Diagonal(ones(5)))
function Gz_func(q)
	cat(G_func(q), tmp, dims=(1, 2))
end

r = residual(z, θ, κ)
∇r = Symbolics.jacobian(r, z) * Gz_func(z)

r_func = eval(Symbolics.build_function(r, z, θ, κ)[2])
∇r_func = eval(Symbolics.build_function(∇r, z, θ)[2])

rz = similar(∇r, Float64)

# options
opts = ContactControl.InteriorPointOptions(diff_sol = false)

∇r_func(rz, ones(nz), ones(nθ))

# Rn + quaternion space
rq_space = rn_quaternion_space(11, x -> Gz_func(x),
			collect([(1:3)..., (8:12)...]),
			collect([(1:3)..., (7:11)...]),
			[collect(4:7)],
			[collect(4:6)])

p0 = [0.0; 0.5; 0.0]
_quat0 = UnitQuaternion(AngleAxis(-0.5 * π, 1.0, 0.0, 0.0))
quat0 = [Rotations.scalar(_quat0); Rotations.vector(_quat0)...]
q0 = [p0; quat0]

p1 = [0.0; 0.5; 0.0]
quat1 = copy(quat0)
q1 = [p1; quat1]

h = 0.1
θ0 = [q0; q1; 0.1]
z0 = copy([q1; 0.001 * randn(5)])

# solver
ip = ContactControl.interior_point(z0, θ0,
	s = rq_space,
	idx_ineq = collect(1:0),
	r! = r_func, rz! = ∇r_func,
	rz = rz,
	opts = opts)

# solve
T = 10
q_hist = [q0, q1]

for t = 1:T-2
	ip.θ .= [q_hist[end-1]; q_hist[end]; h]
	ip.z .= copy([q_hist[end]; 0.001 * randn(5)])
	status = ContactControl.interior_point_solve!(ip)
	push!(q_hist, ip.z[1:nq])
end

function visualize!(vis, model::Body, q;
        Δt = 0.1)

	default_background!(vis)

	r = 0.5

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
visualize!(vis, link, q_hist, Δt = h)
