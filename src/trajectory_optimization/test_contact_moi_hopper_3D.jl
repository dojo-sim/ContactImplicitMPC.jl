include(joinpath(pwd(), "src/trajectory_optimization/moi.jl"))
include(joinpath(pwd(), "src/trajectory_optimization/utils.jl"))

s = get_simulation("hopper_3D", "flat_3D_nc", "flat_nc")

n = s.model.dim.q
m = s.model.dim.u

T = 10
h = 0.1

q0 = [0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.5]
q1 = [0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.5]
qT = [0.25; 0.25; 0.5; 0.0; 0.0; 0.0; 0.5]
q_ref = [0.25; 0.25; 0.5; 0.0; 0.0; 0.0; 0.25]

nz = n * T + m * (T - 2)
np = n * (T - 2)

qq = linear_interpolation(q0, qT, T)
u0 = [0.1 * randn(m) for t = 1:T-2]
z0 = vcat(qq[1], [[qq[t+1]; u0[t]] for t = 1:T-2]..., qq[T])

slack = 1.0

struct Dynamics{T}
	s::Simulation
	ip_dyn::InteriorPoint
	ip_jac::InteriorPoint
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

	# # Rn + quaternion space
	# rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	# 			collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	# 			collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	# 			collect((4:7)),
	# 			collect((4:6)))

	ip_dyn = interior_point(z, θ,
		# s = rq_space,
		idx_ineq = inequality_indices(s.model, s.env),
		# idx_soc = index_soc(s.model, s.env),
		r! = s.res.r!,
		rz! = s.res.rz!,
		rθ! = s.res.rθ!,
		rz = s.rz,
		rθ = s.rθ,
		opts = ip_opts)

	ip_dyn.opts.diff_sol = false

	ip_jac = interior_point(z, θ,
		# s = rq_space,
		idx_ineq = inequality_indices(s.model, s.env),
		# idx_soc = index_soc(s.model, s.env),
		r! = s.res.r!,
		rz! = s.res.rz!,
		rθ! = s.res.rθ!,
		rz = s.rz,
		rθ = s.rθ,
		opts = ip_opts)

	ip_jac.opts.diff_sol = true

	Dynamics(s, ip_dyn, ip_jac, h)
end

d = gen_dynamics(s, h, ip_opts = InteriorPointOptions{Float64}(κ_tol = 0.1 * slack, κ_init = slack))

function f!(d::Dynamics, q0, q1, u1, mode = :dynamics)
	s = d.s
	ip = (mode == :dynamics ? d.ip_dyn : d.ip_jac)
	h = d.h

	z_initialize!(ip.z, s.model, s.env, copy(q1))
	θ_initialize!(ip.θ, s.model, copy(q0), copy(q1), copy(u1), zeros(s.model.dim.w), s.model.μ_world, h)

	# ip.opts.diff_sol = true
	status = interior_point_solve!(ip)

	!status && (@warn "dynamics failure")
end

function f(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1, :dynamics)
	return copy(d.ip_dyn.z[1:d.s.model.dim.q])
end

f(d, q0, q1, zeros(m))

function fq0(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1, :jacobian)
	return copy(d.ip_jac.δz[1:d.s.model.dim.q, 1:d.s.model.dim.q])
end

fq0(d, q0, q1, zeros(m))

function fq1(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1, :jacobian)
	return copy(d.ip_jac.δz[1:d.s.model.dim.q, d.s.model.dim.q .+ (1:d.s.model.dim.q)])
end

fq1(d, q0, q1, zeros(m))

function fu1(d::Dynamics, q0, q1, u1)
	f!(d, q0, q1, u1, :jacobian)
	return copy(d.ip_jac.δz[1:d.s.model.dim.q, 2 * d.s.model.dim.q .+ (1:d.s.model.dim.u)])
end

fu1(d, q0, q1, zeros(m))

function obj_config(q)
	Qt = Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 100.0])
	return transpose(q - q_ref) * Qt * (q - q_ref)
end

function obj_configT(q)
	QT = Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 100.0])
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
	Rt = Diagonal(0.001 * ones(m))
	return transpose(ut) * Rt * ut
end

@variables u_sym[1:m]

obj_ctrl_sym = obj_ctrl(u_sym)
obj_ctrl_sym = simplify.(obj_ctrl_sym)

obj_ctrl_u_sym = Symbolics.gradient(obj_ctrl_sym, u_sym, simplify = true)

obj_ctrl_func = eval(Symbolics.build_function(obj_ctrl_sym, u_sym))
obj_ctrl_u_func = eval(Symbolics.build_function(obj_ctrl_u_sym, u_sym)[1])

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
		J += obj_ctrl_func(ut)
	end

	return J
end

z0
moi_obj(z0)

function ∇moi_obj!(g, z)
	for t = 1:T
		qt = z[q_idx[t]]
		g[q_idx[t]] .= t == T ? obj_configT_func(qt) : obj_config_q_func(qt)
	end

	for t = 1:T-2
		ut = z[u_idx[t]]
		g[u_idx[t]] .= obj_ctrl_u_func(ut)
	end

	nothing
end

g0 = zeros(nz)
∇moi_obj!(g0, z0)
g0

function moi_con!(c, z)
	for t = 1:T-2
		q0 = z[q_idx[t]]
		q1 = z[q_idx[t+1]]
		q2 = z[q_idx[t+2]]
		u1 = z[u_idx[t]]
		c[(t - 1) * n .+ (1:n)] .= f(d, q0, q1, u1) - q2
	end
	nothing
end

c0 = zeros(np)
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
	nothing
end

j0 = zeros((n * n + n * n + n * m + n * n) * (T - 2))
∇moi_con!(j0, z0)

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

    return collect(zip(row, col))
end

sparsity_jacobian(nz, np)

zl, zu = primal_bounds(nz)
zl[1:n] = q0
zu[1:n] = q0
zl[n .+ (1:n)] = q1
zu[n .+ (1:n)] = q1

zl[end-n+1:end] = qT
zu[end-n+1:end] = qT

for t = 1:T-2
	zl[u_idx[t]] = [-10.0; -10.0; -10.0]
	zu[u_idx[t]] = [10.0; 10.0; 10.0]
end

for t = 3:T-1
	zl[q_idx[t][7]] = 0.0
	zu[q_idx[t][7]] = 1.0
end

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
slack = 0.1 * slack
cl, cu = constraint_bounds(np)
cl .= -slack
cu .= slack
d.ip_dyn.opts.κ_init = slack
d.ip_dyn.opts.κ_tol = 0.1 * slack
d.ip_jac.opts.κ_init = slack
d.ip_jac.opts.κ_tol = 0.1 * slack
prob = ProblemMOI(nz, np,
    sparsity_jac=sparsity_jacobian(nz, np),
    primal_bounds=(zl, zu),
    constraint_bounds=(cl, cu),
    hessian_lagrangian=false)

z_sol = solve(z_sol + 0.001 * randn(nz), prob,
        nlp=:ipopt,
        tol=1.0e-3,
        c_tol=1.0e-3,
        max_iter=1000)

c0 = zeros(np)
moi_con!(c0, z_sol)
norm(c0, Inf)
@show d.ip_dyn.κ
@show d.ip_jac.κ

x_traj = [z_sol[q_idx[t]] for t = 1:T]
u_traj = [z_sol[u_idx[t]] for t = 1:T-2]

plot(hcat(x_traj...)')
plot(hcat(u_traj...)', linetype = :steppost)

include(joinpath(@__DIR__, "..", "dynamics", "hopper_3D", "visuals.jl"))
vis = Visualizer()
render(vis)
anim = visualize_robot!(vis, s.model, x_traj, h = h)
