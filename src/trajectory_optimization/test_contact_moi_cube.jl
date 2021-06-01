include(joinpath(pwd(), "src/trajectory_optimization/moi.jl"))
include(joinpath(pwd(), "src/trajectory_optimization/utils.jl"))
include(joinpath(pwd(), "src/dynamics/box/model.jl"))

s = get_simulation("box", "flat_3D_nc", "flat_nc")
s.model.μ_world = 0.5

n = s.model.dim.q
m = s.model.dim.u

T = 10
h = 0.1

q0 = [0.5; 0.5; 0.5; 1.0; 0.0; 0.0; 0.0]
q1 = [0.5; 0.5; 0.5; 1.0; 0.0; 0.0; 0.0]
quat1 = UnitQuaternion(RotY(0.0) * RotX(pi / 4.0))

# quat1 = UnitQuaternion(RotY(-1.0 * atan(1.0 / sqrt(2.0))) * RotX(pi / 4.0))
qT = [0.5; 0.0; 0.5 * sqrt(2.0); quat1.w; quat1.x; quat1.y; quat1.z]
q_ref = copy(qT)

k = kinematics(s.model, qT)
p4 = k[3 * 3 .+ (1:3)]
p8 = k[7 * 3 .+ (1:3)]

nz = n * T + m * (T - 2)
np = n * (T - 2)# + 2 * 3 * T

z0 = rand(nz)

slack = 0.1

struct Dynamics{T}
	s::Simulation
	ip::InteriorPoint
	h::T
end

function gen_dynamics(s::Simulation, h;
		ip_opts =  InteriorPointOptions{Float64}(
						r_tol = 1.0e-5,
						κ_tol = 0.1,
						κ_init = 0.1,
						diff_sol = true))

	z = zeros(num_var(s.model, s.env))
	θ = zeros(num_data(s.model))

	# Rn + quaternion space
	rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
				collect([(1:3)..., (8:num_var(s.model, s.env))...]),
				collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
				collect((4:7)),
				collect((4:6)))

	ip = interior_point(z, θ,
		s = rq_space,
		idx_ineq = inequality_indices(s.model, s.env),
		idx_soc = soc_indices(s.model, s.env),
		r! = s.res.r!,
		rz! = s.res.rz!,
		rθ! = s.res.rθ!,
		rz = s.rz,
		rθ = s.rθ,
		opts = ip_opts)

	Dynamics(s, ip, h)
end

d = gen_dynamics(s, h, ip_opts = InteriorPointOptions{Float64}(κ_tol = slack, κ_init = slack))

function f!(d::Dynamics, q0, q1, u1)
	s = d.s
	ip = d.ip
	h = d.h

	z_initialize!(ip.z, s.model, s.env, copy(q1))
	θ_initialize!(ip.θ, s.model, copy(q0), copy(q1), copy(u1), zeros(s.model.dim.w), s.model.μ_world, h)

	ip.opts.diff_sol = true
	status = interior_point!(ip)

	!status && (@warn "dynamics failure")
end

function f(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1)
	return copy(d.ip.z[1:d.s.model.dim.q])
end

f(d, q0, q1, zeros(m))

function fq0(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1)
	return copy(d.ip.δz[1:d.s.model.dim.q, 1:d.s.model.dim.q])
end

fq0(d, q0, q1, zeros(m))

function fq1(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1)
	return copy(d.ip.δz[1:d.s.model.dim.q, d.s.model.dim.q .+ (1:d.s.model.dim.q)])
end

fq1(d, q0, q1, zeros(m))

function fu1(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1)
	return copy(d.ip.δz[1:d.s.model.dim.q, 2 * d.s.model.dim.q .+ (1:d.s.model.dim.u)])
end

fu1(d, q0, q1, zeros(m))

function obj_config(q)
	Qt = Diagonal(1.0e-1 * [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0])
	return transpose(q - q_ref) * Qt * (q - q_ref)
end

function obj_configT(q)
	QT = Diagonal(1.0e-1 * [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0])
	return transpose(q - q_ref) * QT * (q - q_ref)
end

@variables q_sym[1:n]

obj_config_sym = obj_config(q_sym)
obj_config_sym = simplify.(obj_config_sym)

obj_config_q_sym = Symbolics.gradient(obj_config_sym, q_sym, simplify = true)

obj_config_func = eval(Symbolics.build_function(obj_config_sym, q_sym))
obj_config_q_func = eval(Symbolics.build_function(obj_config_q_sym, q_sym)[1])

obj_configT_sym = obj_configT(q_sym)
obj_configT_sym = simplify.(obj_configT_sym)

obj_configT_q_sym = Symbolics.gradient(obj_configT_sym, q_sym, simplify = true)

obj_configT_func = eval(Symbolics.build_function(obj_configT_sym, q_sym))
obj_configT_q_func = eval(Symbolics.build_function(obj_configT_q_sym, q_sym)[1])

function obj_ctrl(ut)
	Rt = Diagonal(0.1 * ones(m))
	return transpose(ut) * Rt * ut
end

@variables u_sym[1:m]

obj_ctrl_sym = obj_ctrl(u_sym)
obj_ctrl_sym = simplify.(obj_ctrl_sym)

obj_ctrl_u_sym = Symbolics.gradient(obj_ctrl_sym, u_sym, simplify = true)

obj_ctrl_func = eval(Symbolics.build_function(obj_ctrl_sym, u_sym))
obj_ctrl_u_func = eval(Symbolics.build_function(obj_ctrl_u_sym, u_sym)[1])

@variables q0_sym[1:n]
@variables h_sym[1:1]

function obj_vel(q0, q1, h)
	Qt = Diagonal(0.001 * ones(n))
	return transpose(q1 - q0) * Qt * (q1 - q0) / h[1]^2.0
end

obj_vel_sym = obj_vel(q0_sym, q_sym, h_sym)
obj_vel_sym = simplify.(obj_vel_sym)

obj_vel_q0_sym = Symbolics.gradient(obj_vel_sym, q0_sym, simplify = true)
obj_vel_q1_sym = Symbolics.gradient(obj_vel_sym, q_sym, simplify = true)

obj_vel_func = eval(Symbolics.build_function(obj_vel_sym, q0_sym, q_sym, h_sym))
obj_vel_q0_func = eval(Symbolics.build_function(obj_vel_q0_sym, q0_sym, q_sym, h_sym)[1])
obj_vel_q1_func = eval(Symbolics.build_function(obj_vel_q1_sym, q0_sym, q_sym, h_sym)[1])

q_idx = [(1:n), [n + (t - 1) * (n + m) .+ (1:n) for t = 1:T-2]..., n + (T - 1 - 1) * (n + m) .+ (1:n)]
u_idx = [n + (t - 1) * (n + m) + n .+ (1:m) for t = 1:T-2]

function moi_obj(z)
	J = 0.0

	for t = 1:T
		qt = z[q_idx[t]]
		J += t == T ? obj_configT_func(qt) : obj_config_func(qt)
	end

	for t = 1:T-2
		ut = z[u_idx[t]]
		J += obj_config_func(ut)
	end

	for t = 1:T-1
		q0 = z[q_idx[t]]
		q1 = z[q_idx[t+1]]

		J += obj_vel_func(q0, q1, h)
	end

	return J
end

obj_vel_func(rand(7), rand(7), h)

z0
moi_obj(z0)

function ∇moi_obj!(g, z)
	g .= 0.0

	for t = 1:T
		qt = z[q_idx[t]]
		g[q_idx[t]] .+= t == T ? obj_configT_func(qt) : obj_config_q_func(qt)
	end

	for t = 1:T-2
		ut = z[u_idx[t]]
		g[u_idx[t]] .+= obj_ctrl_u_func(ut)
	end

	for t = 1:T-1
		q0 = z[q_idx[t]]
		q1 = z[q_idx[t+1]]

		g[q_idx[t]] .+= obj_vel_q0_func(q0, q1, h)
		g[q_idx[t+1]] .+= obj_vel_q1_func(q0, q1, h)
	end

	nothing
end

g0 = zeros(nz)
∇moi_obj!(g0, z0)
g0

function moi_con!(c, z)
	c .= 0.0
	for t = 1:T-2
		q0 = z[q_idx[t]]
		q1 = z[q_idx[t+1]]
		q2 = z[q_idx[t+2]]
		u1 = z[u_idx[t]]
		c[(t - 1) * n .+ (1:n)] .= f(d, q0, q1, u1) - q2
	end

	# for t = 1:T
	# 	qt = z[q_idx[t]]
	# 	k = kinematics(model, qt)
	# 	c[n * (T - 2) + (t - 1) * 2 * 3 .+ (1:3)] .= k[3 * 3 .+ (1:3)] - p4
	# 	c[n * (T - 2) + (t - 1) * 2 * 3 + 3 .+ (1:3)] .= k[7 * 3 .+ (1:3)] - p8
	# end

	nothing
end

qq = linear_interpolation(q0, qT, T)
u0 = [0.001 * randn(m) for t = 1:T-2]
z0 = vcat(qq[1], [[qq[t+1]; u0[t]] for t = 1:T-2]..., qq[T])


c0 = zeros(np)
z0[4:7] = [1.0; 0.0; 0.0; 0.0]
moi_con!(c0, z0)
c0

function ∇moi_con!(j, z)
	shift = 0
	for t = 1:T-2
		q0 = z[q_idx[t]]
		q1 = z[q_idx[t+1]]
		q2 = z[q_idx[t+2]]
		u1 = z[u_idx[t]]

		r_idx = (t - 1) * n .+ (1:n)

		c_idx = q_idx[t]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(fq0(d, q0, q1, u1))
		shift += s

		c_idx = q_idx[t+1]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(fq1(d, q0, q1, u1))
		shift += s

		c_idx = u_idx[t]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(fu1(d, q0, q1, u1))
		shift += s

		c_idx = q_idx[t+2]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(Diagonal(-1.0 * ones(n)))
		shift += s
	end

	# for t = 1:T
	# 	q = z[q_idx[t]]
	# 	k = kinematics(model, q)
	# 	k_func(w) = kinematics(model, w)
	#
	# 	# c[n * (T - 2) + (t - 1) * 2 * 3 .+ (1:3)] = k[3 * 3 .+ (1:3)] - p4
	# 	r_idx = collect(n * (T - 2) + (t - 1) * 2 * 3 .+ (1:3))
	# 	c_idx = q_idx[t]
	# 	s = length(r_idx) * length(c_idx)
	# 	k4(w) = kinematics(model, w)[3 * 3 .+ (1:3)]
	# 	j[shift .+ (1:s)] .= vec(ForwardDiff.jacobian(k4, q))
	#
	# 	# c[n * (T - 2) + (t - 1) * 2 * 3 + 3 .+ (1:3)] = k[6 * 3 .+ (1:3)] - p7
	# 	r_idx = collect(n * (T - 2) + (t - 1) * 2 * 3 + 3 .+ (1:3))
	# 	c_idx = q_idx[t]
	# 	s = length(r_idx) * length(c_idx)
	# 	k8(w) = kinematics(model, w)[7 * 3 .+ (1:3)]
	# 	j[shift .+ (1:s)] .= vec(ForwardDiff.jacobian(k8, q))
	# end

	nothing
end

j0 = zeros((n * n + n * n + n * m + n * n) * (T - 2) + 0.0 * (2 * 3 * n) * T)
∇moi_con!(j0, z0)

j0

# user can overwrite sparsity_jacobian and sparsity_hessian
function sparsity_jacobian(nz, np)

    row = []
    col = []

	for t = 1:T-2
		r_idx = (t - 1) * n .+ (1:n)

		c_idx = q_idx[t]
		row_col!(row,col,r_idx,c_idx)

		c_idx = q_idx[t+1]
		row_col!(row,col,r_idx,c_idx)

		c_idx = u_idx[t]
		row_col!(row,col,r_idx,c_idx)

		c_idx = q_idx[t+2]
		row_col!(row,col,r_idx,c_idx)
	end

	# for t = 1:T
	# 	# q = z[q_idx[t]]
	# 	# k = kinematics(model, q)
	# 	# k_func(w) = kinematics(model, w)
	#
	# 	# c[n * (T - 2) + (t - 1) * 2 * 3 .+ (1:3)] = k[3 * 3 .+ (1:3)] - p4
	# 	r_idx = n * (T - 2) + (t - 1) * 2 * 3 .+ (1:3)
	# 	c_idx = q_idx[t]
	# 	row_col!(row,col,r_idx,c_idx)
	#
	# 	# c[n * (T - 2) + (t - 1) * 2 * 3 + 3 .+ (1:3)] = k[6 * 3 .+ (1:3)] - p7
	# 	r_idx = n * (T - 2) + (t - 1) * 2 * 3 + 3 .+ (1:3)
	# 	c_idx = q_idx[t]
	# 	row_col!(row,col,r_idx,c_idx)
	#
	# end

    return collect(zip(row, col))
end

sparsity_jacobian(nz, np)

zl, zu = primal_bounds(nz)
zl[1:n] = q0 .- slack
zu[1:n] = q0 .+ slack
zl[n .+ (1:n)] = q1 .- slack
zu[n .+ (1:n)] = q1 .+ slack

zl[end-n+1:end] = qT .- slack
zu[end-n+1:end] = qT .+ slack

# for t = 1:T-2
# 	zl[u_idx[t]] = [-10.0; -10.0; -10.0]
# 	zu[u_idx[t]] = [10.0; 10.0; 10.0]
# end

cl, cu = constraint_bounds(np)
cl .= -slack
cu .= slack

# j0 = zeros((n * n + n * n + n * m + n * n) * (T - 2))
# ∇moi_con!(j0, z0)
# J0 = zeros(np, nz)
# spar = sparsity_jacobian(nz, np)
# for (i, idx) in enumerate(spar)
# 	J0[idx...] = j0[i]
# end

# show(J0)

# plot(hcat(qq...)')
prob = ProblemMOI(nz, np,
    sparsity_jac=sparsity_jacobian(nz, np),
    primal_bounds=(zl, zu),
    constraint_bounds=(cl, cu),
    hessian_lagrangian=false)

z_sol = solve(z0, prob,
        nlp=:ipopt,
        tol=1.0e-2,
        c_tol=1.0e-2,
        max_iter=1000)

c0 = zeros(np)
moi_con!(c0, z_sol)
norm(c0)

x_traj = [z_sol[q_idx[t]] for t = 1:T]
u_traj = [z_sol[u_idx[t]] for t = 1:T-2]

plot(hcat(x_traj...)')
plot(hcat(u_traj...)', linetype = :steppost)

##
slack = 0.5 * slack

zl, zu = primal_bounds(nz)
zl[1:n] = q0 .- slack
zu[1:n] = q0 .+ slack
zl[n .+ (1:n)] = q1 .- slack
zu[n .+ (1:n)] = q1 .+ slack

zl[end-n+1:end] = qT .- slack
zu[end-n+1:end] = qT .+ slack


# for t = 1:T-2
# 	zl[u_idx[t]] = [-10.0; -10.0; -10.0]
# 	zu[u_idx[t]] = [10.0; 10.0; 10.0]
# end

cl, cu = constraint_bounds(np)
cl .= -slack
cu .= slack
d.ip.opts.κ_init = slack
d.ip.opts.κ_tol = slack
prob = ProblemMOI(nz, np,
    sparsity_jac=sparsity_jacobian(nz, np),
    primal_bounds=(zl, zu),
    constraint_bounds=(cl, cu),
    hessian_lagrangian=false)

z_sol = solve(z_sol + 0.001 * randn(nz), prob,
        nlp=:ipopt,
        tol=1.0e-2,
        c_tol=1.0e-2,
        max_iter=1000)

c0 = zeros(np)
moi_con!(c0, z_sol)
norm(c0, Inf)
@show d.ip.κ

x_traj = [z_sol[q_idx[t]] for t = 1:T]
u_traj = [z_sol[u_idx[t]] for t = 1:T-2]

plot(hcat(x_traj...)')
plot(hcat(u_traj...)', linetype = :steppost)

include(joinpath(@__DIR__, "..", "dynamics", "box", "visuals.jl"))
vis = Visualizer()
render(vis)
anim = visualize!(vis, s.model, x_traj, Δt = h)
