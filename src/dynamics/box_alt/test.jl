include(joinpath(@__DIR__, "..", "..", "dynamics/quaternions.jl"))

function visualize!(vis, p, q; Δt = 0.1)
	setvisible!(vis["/Background"], true)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	setvisible!(vis["/Grid"], false)

    setobject!(vis[:satellite],
    	Rect(Vec(-0.25, -0.25, -0.25),Vec(0.5, 0.5, 0.5)),
    	MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    arrow_x = ArrowVisualizer(vis[:satellite][:arrow_x])
    mat = MeshPhongMaterial(color=RGBA(1.0, 0.0, 0.0, 1.0))
    setobject!(arrow_x, mat)
    settransform!(arrow_x,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.75, 0.0, 0.0),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    arrow_y = ArrowVisualizer(vis[:satellite][:arrow_y])
    mat = MeshPhongMaterial(color=RGBA(0.0, 1.0, 0.0, 1.0))
    setobject!(arrow_y, mat)
    settransform!(arrow_y,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.0, 0.75, 0.0),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    arrow_z = ArrowVisualizer(vis[:satellite][:arrow_z])
    mat = MeshPhongMaterial(color=RGBA(0.0, 0.0, 1.0, 1.0))
    setobject!(arrow_z, mat)
    settransform!(arrow_z,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.0, 0.0, 0.75),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

     # for t = 1:length(q)
    	#  MeshCat.atframe(anim, t) do
    	# 	 settransform!(vis["satellite"],
    	# 		   compose(Translation((q[t][1:3] + [-0.25; -0.25; -0.25])...),
    	# 				 LinearMap(quaternion_rotation_matrix(q[t][4:7]))))
    	#  end
     # end

	 for t = 1:length(q)
		MeshCat.atframe(anim, t) do
			settransform!(vis["satellite"],
				  compose(Translation(([-0.0; -0.0; -0.0])...),
						LinearMap(UnitQuaternion(q[t][1:4]))))
		end
	 end
	#
	#
    MeshCat.setanimation!(vis, anim)
end

vis = Visualizer()
# render(vis)
open(vis)

nq = 4
nv = 3
H = [zeros(1, 3); I]

J = Diagonal([1.0; 2.0; 3.0])

h = 0.01
T = 1000

function f(x, u)
	q = x[1:4]
	ω = x[5:7]
	τ = u[1:3]

	[0.5 * L_multiply(q) * H * ω;
	 J \ (τ - cross(ω, J * ω))]
 end

function f(x, u, h)
    x + h * f(x + 0.5 * h * f(x, u), u) # explicit midpoint
end

q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 10.0; 0.1]

x1 = [q1; ω1]
u1 = zeros(3)

x = [copy(x1)]
for t = 1:T
	push!(x, f(x[end], u1, h))
end

visualize!(vis, nothing, x, Δt = h)

# implicit
function f(x⁺, x, u, h)
	x⁺ - (x + h * f(0.5 * (x + x⁺), u)) # implciit midpoint
end

function newton(x, u, h)
	x⁺ = copy(x)
	_r(z) = f(z, x, u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		x⁺[1:4] = L_multiply(x[1:4]) * φ(α * Δ[1:3])
		x⁺[5:7] += α * Δ[4:6]

		return x⁺
	end

	r = _r(x⁺)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, x⁺) * [attitude_jacobian(x⁺) zeros(4, 3); zeros(3,3) I]
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		x_cand = cand(x⁺, Δ, α)
		r_cand = _r(x_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) >= norm(r)
			α *= 0.5
			x_cand = cand(x⁺, Δ, α)
			r_cand = _r(x_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		x⁺ = x_cand
		r = r_cand

		if norm(r) < 1.0e-8
			return x⁺
		end
	end
	@warn "newton failure"
	return x
end

newton(x1, u1, h)

x = [copy(x1)]
for t = 1:T
	push!(x, newton(x[end], u1, h))
end

visualize!(vis, nothing, x, Δt = h)

function newton_variational(x, u, h)
	q2 = copy(x)[1:4]
	ω2 = copy(x)[5:7]

	function f(ω1, ω2, q2, u, h)
		(J * ω2 * sqrt(4.0 / h^2.0 - ω2' * ω2)
			+ cross(ω2, J * ω2)
			- J * ω1 * sqrt(4.0 / h^2.0 - ω1' * ω1)
			+ cross(ω1, J * ω1)
			- 2.0 * u[1:3])
	end

	_r(z) = f(x[5:7], z, x[1:4], u, h)

	function cand(x, Δ, α)
		x⁺ = copy(x)
		# x⁺[1:4] = L_multiply(x[1:4]) * φ(α * Δ[1:3])
		x⁺[1:3] += α * Δ[1:3]

		return x⁺
	end

	r = _r(ω2)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, ω2)# * [attitude_jacobian(x⁺) zeros(4, 3); zeros(3,3) I]
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		ω_cand = cand(ω2, Δ, α)
		r_cand = _r(ω_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) >= norm(r)
			α *= 0.5
			x_cand = cand(ω2, Δ, α)
			r_cand = _r(ω_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		ω2 = ω_cand
		r = r_cand

		if norm(r) < 1.0e-10
			q3 = 0.5 * h * L_multiply(q2) * [sqrt((2.0 / h)^2.0 - ω2' * ω2); ω2]
			return [q3; ω2]
		end
	end
	@warn "newton failure"
	return x
end

newton_variational(x1, u1, h)

x = [copy(x1)]
for t = 1:T
	push!(x, newton_variational(x[end], u1, h))
end

visualize!(vis, nothing, x, Δt = h)

q_variational = [xt[1:4] for xt in x]
w_variational = [xt[5:7] for xt in x]
plot(hcat(w_variational...)')
plot(hcat(q_variational...)')
plot([norm(q) for q in q_variational])

function newton_variational_2(q1, q2, u, h)
	q3 = copy(q2)

	# @show ω_finite_difference(q1, q2, h)
	# @show norm(q3)
	# @show q3
	function f(q1, q2, q3, u, h)
		ω1 = ω_finite_difference(q1, q2, h)
		ω2 = ω_finite_difference(q2, q3, h)

		(J * ω2 * sqrt(4.0 / h^2.0 - ω2' * ω2)
			+ cross(ω2, J * ω2)
			- J * ω1 * sqrt(4.0 / h^2.0 - ω1' * ω1)
			+ cross(ω1, J * ω1)
			- 2.0 * u[1:3])
	end

	_r(z) = f(q1, q2, z, u, h)

	function cand(q3, Δ, α)
		return L_multiply(q3) * φ(α * Δ[1:3])
	end

	r = _r(q3)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, q3) * attitude_jacobian(q3)
		Δ = -(∇r' * ∇r) \ (∇r' * r)
		α = 1.0

		q3_cand = cand(q3, Δ, α)
		r_cand = _r(q3_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) > norm(r)
			α *= 0.5
			q3_cand = cand(q3, Δ, α)
			r_cand = _r(q3_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		q3 = q3_cand
		r = r_cand

		# @show ω_finite_difference(q2, q3, h)
		# @show norm(q3)
		# @show q3
		if norm(r) < 1.0e-10
			return q3
		end
	end
	@warn "newton failure"
	return q2
end

x = [copy(q1), copy(Array(0.5 * h * L_multiply(q1) * [sqrt((2.0 / h)^2.0 - ω1' * ω1); ω1]))]

for t = 1:500
	push!(x, newton_variational_2(x[end-1], x[end], u1, h))
end

visualize!(vis, nothing, x, Δt = h)

function newton_variational_3(q1, q2, u, h)
	q3 = copy(q2)

	# @show ω_finite_difference(q1, q2, h)
	# @show norm(q3)
	# @show q3
	function f(q1, q2, q3, u, h)
		v1 = (q2[1:3] - q1[1:3]) ./ h
		v2 = (q3[1:3] - q2[1:3]) ./ h

		ω1 = ω_finite_difference(q1[4:7], q2[4:7], h)
		ω2 = ω_finite_difference(q2[4:7], q3[4:7], h)

		D1L1, D2L1 = [0.0; 0.0; -9.81], v1
		D1L2, D2L2 = [0.0; 0.0; -9.81], v2

		[0.5 * h * D1L1 + D2L1 + 0.5 * h * D1L2 - D2L2;
			(J * ω2 * sqrt(4.0 / h^2.0 - ω2' * ω2)
			+ cross(ω2, J * ω2)
			- J * ω1 * sqrt(4.0 / h^2.0 - ω1' * ω1)
			+ cross(ω1, J * ω1)
			- 2.0 * u[1:3])]
	end

	_r(z) = f(q1, q2, z, u, h)

	function cand(q3, Δ, α)
		return [q3[1:3] + α .* Δ[1:3];
			    L_multiply(q3[4:7]) * φ(α * Δ[4:6])]
	end

	r = _r(q3)
	@show norm(r)

	for i = 1:100
		∇r = ForwardDiff.jacobian(_r, q3) * [I zeros(3, 3);
		                                      zeros(4, 3) attitude_jacobian(q3[4:7])]
		# Δ = -(∇r' * ∇r) \ (∇r' * r)
		Δ = -∇r \ r

		α = 1.0

		q3_cand = cand(q3, Δ, α)
		r_cand = _r(q3_cand)
		@show norm(r_cand)
		iter = 0

		while norm(r_cand) > norm(r)
			α *= 0.5
			q3_cand = cand(q3, Δ, α)
			r_cand = _r(q3_cand)
			@show norm(r_cand)

			iter += 1
			if iter > 100
				@error "line search failure"
				break
			end
		end

		q3 = q3_cand
		r = r_cand

		# @show ω_finite_difference(q2, q3, h)
		# @show norm(q3)
		# @show q3
		if norm(r) < 1.0e-10
			return q3
		end
	end
	@warn "newton failure"
	return q2
end

p1 = [0.0; 0.0; 0.0]
v1 = [10.0; 0.0; 5.0]
x = [[copy(p1); copy(q1)],
	 [copy(p1) + h * copy(v1); copy(Array(0.5 * h * L_multiply(q1) * [sqrt((2.0 / h)^2.0 - ω1' * ω1); ω1]))]]

for t = 1:100
	push!(x, newton_variational_3(x[end-1], x[end], u1, h))
end
include(joinpath(module_dir(), "src/dynamics/rigid_body/visuals.jl"))

visualize!(vis, nothing, x, Δt = h)


struct RQuat <: Space
	n::Int
	r_idx
	Δr_idx
	quat_idx
	Δquat_idx
end

function rquat(dim, r_idx, Δr_idx, quat_idx, Δquat_idx)
	RQuat(dim, r_idx, Δr_idx, quat_idx, Δquat_idx)
end

rquat_space = rquat(6,
			collect((1:3)),
			collect((1:3)),
			collect((4:7)),
			collect((4:6)))

function candidate_point!(z̄::Vector{T}, s::RQuat, z::Vector{T}, Δ::Vector{T}, α::T) where T
    z̄[s.r_idx] .= z[s.r_idx] - α .* Δ[s.Δr_idx]
	z̄[s.quat_idx] .= L_multiply(z[s.quat_idx]) * φ(-1.0 * α .* Δ[s.Δquat_idx])
	return nothing
end

function f(q1, q2, q3, u, h)
	v1 = (q2[1:3] - q1[1:3]) ./ h
	v2 = (q3[1:3] - q2[1:3]) ./ h

	ω1 = ω_finite_difference(q1[4:7], q2[4:7], h)
	ω2 = ω_finite_difference(q2[4:7], q3[4:7], h)

	D1L1, D2L1 = [0.0; 0.0; -9.81], v1
	D1L2, D2L2 = [0.0; 0.0; -9.81], v2

	[0.5 * h * D1L1 + D2L1 + 0.5 * h * D1L2 - D2L2;
		-1.0 * (J * ω2 * sqrt(4.0 / h^2.0 - ω2' * ω2)
		+ cross(ω2, J * ω2)
		- J * ω1 * sqrt(4.0 / h^2.0 - ω1' * ω1)
		+ cross(ω1, J * ω1)
		- 2.0 * u[1:3])]
end

function G_func(x)
	q = x[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(q)]
end

function r!(r, z, θ, κ)
	q1 = θ[1:7]
	q2 = θ[7 .+ (1:7)]
	u1 = θ[2 * 7 .+ (1:3)]
	h = θ[2 * 7 + 3 + 1]
	r .= f(q1, q2, z, u1, h)
	nothing
end

function rz!(rz, z, θ)
	q1 = θ[1:7]
	q2 = θ[7 .+ (1:7)]
	u1 = θ[2 * 7 .+ (1:3)]
	h = θ[2 * 7 + 3 + 1]
	_f(y) = f(q1, q2, y, u1, h)
	rz .= ForwardDiff.jacobian(_f, z) * G_func(z)
end

rzz = zeros(6,6)
rz!(rzz, rand(7), rand(18))

# options
opts = ContactControl.InteriorPointOptions(
	r_tol = 1.0e-6, κ_tol = 1.0e-6,
	diff_sol = false,
	solver = :lu_solver)
#
# rr = zeros(6)
# r!(rr, z, θ, 0.0)
# rr
p1 = [0.0; 0.0; 0.0]
v1 = [10.0; 0.0; 5.0]
ω1 = [0.0; 10.0; 0.1]
x = [[copy(p1); copy(q1)],
	 [copy(p1) + h * copy(v1); copy(Array(0.5 * h * L_multiply(q1) * [sqrt((2.0 / h)^2.0 - ω1' * ω1); ω1]))]]

z = copy(x[end])
θ = [x[end-1]; x[end]; u1; h]

# solver
ip = ContactControl.interior_point(z, θ,
	s = rquat_space,
	r! = r!, rz! = rz!,
	rz = zeros(6, 6),
	opts = opts)

# solve
for t = 1:T
	ip.z .= copy(x[end])
	ip.θ .= [x[end-1]; x[end]; u1; h]

	status = ContactControl.interior_point!(ip)

	push!(x, copy(ip.z))
end

include(joinpath(module_dir(), "src/dynamics/rigid_body/visuals.jl"))

visualize!(vis, nothing, x, Δt = h)
