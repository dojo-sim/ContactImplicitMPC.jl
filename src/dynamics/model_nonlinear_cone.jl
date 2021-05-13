function unpack_z(model::ContactDynamicsModel, z)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b

	# system variables
	off = 0
	q2 =  z[off .+ (1:nq)]
	off += nq
	γ1 =  z[off .+ (1:nc)]
	off += nc
	b1 =  z[off .+ (1:nb)]
	off += nb
	# ψ1 =  z[off .+ (1:nc)]
	# off += nc
	η1 =  z[off .+ (1:(nb + nc))]
	off += (nb + nc)
	s1 = z[off .+ (1:nc)]
	off += nc
	s2 = z[off .+ (1:nc)]
	off += nc
	return q2, γ1, b1, η1, s1, s2
end

function pack_z(model::ContactDynamicsModel, q2, γ1, b1, η1)
	s1 = ϕ_func(model, q2)
	s2 = model.μ_world .* γ1
	return [q2; γ1; b1; η1; s1; s2]
end

function z_initialize!(z, model::ContactDynamicsModel, q1)
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b
	ne = dim(model.env)

    z .= 0.1
    z[1:nq] = q1

	# second-order cone initializations
	z[nq + nc + nb .+ (1:(nb + nc))] = vcat([[1.0; 0.1 * ones(ne - 1)] for i = 1:model.dim.c]...)
	z[nq + nc + nb + nb + nc + nc .+ (1:nc)] .= 1.0

	return z
end

function num_var(model::ContactDynamicsModel)
	dim = model.dim
	dim.q + dim.c + dim.b + (dim.b + dim.c) + 2 * dim.c
end

function num_data(model::ContactDynamicsModel)
	dim = model.dim
	dim.q + dim.q + dim.u + dim.w + 1 + 1
end

function inequality_indices(model::ContactDynamicsModel)
	collect([(model.dim.q .+ (1:model.dim.c))...,
	         (model.dim.q + model.dim.c + model.dim.b + model.dim.b + model.dim.c .+ (1:model.dim.c))...])
end

function soc_indices(model::ContactDynamicsModel)
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b
	ne = dim(model.env)

	b_idx = nq + nc .+ (1:nb)
	η_idx = nq + nc + nb .+ (1:(nb + nc))
	s2_idx = nq + nc + nb + nb + nc + nc .+ (1:nc)

	pr_idx = [[s2_idx[i]; b_idx[(i - 1) * nb .+ (1:nb)]] for i = 1:nc]
	du_idx = [[η_idx[(i - 1) * (nb + 1) .+ (1:(nb + 1))]...] for i = 1:nc]

	[pr_idx..., du_idx...]
end

function contact_forces(model::ContactDynamicsModel, γ1, b1, q2, k)
	nc = model.dim.c
	nb = model.dim.b
	ne = dim(model.env)
	λ1 = vcat([transpose(rotation(model.env, k[(i-1) * (ne - 1) .+ (1:ne)])) * [b1[(i-1) * (ne-1) .+ (1:(ne-1))]; γ1[i]] for i = 1:nc]...) # TODO: make efficient
end

function velocity_stack(model::ContactDynamicsModel, q1, q2, k, h)
	nc = model.dim.c
	ne = dim(model.env)
	v = J_fast(model, q2) * (q2 - q1) / h[1]
	v_surf = [rotation(model.env, k[(i-1) * ne .+ (1:ne)]) * v[(i-1) * ne .+ (1:ne)] for i = 1:nc]
	vT_stack = vcat([[v_surf[i][1:ne-1]] for i = 1:nc]...)
end

function residual(model::ContactDynamicsModel, z, θ, κ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, η1, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)
	k = kinematics(model, q2)
	λ1 = contact_forces(model, γ1, b1, q2, k)
	vT = velocity_stack(model, q1, q2, k, h)

	ne = dim(model.env)

	[model.dyn.d(h, q0, q1, u1, w1, λ1, q2);
	 s1 - ϕ;
	 vcat([η1[(i - 1) * ne .+ (2:ne)] - vT[(i - 1) * (ne - 1) .+ (1:(ne - 1))] for i = 1:model.dim.c]...);
	 s2 - μ[1] * γ1;
	 γ1 .* s1 .- κ;
	 vcat([second_order_cone_product(η1[(i - 1) * ne .+ (1:ne)], [s2[i]; b1[(i-1) * (ne - 1) .+ (1:(ne - 1))]]) - [κ; zeros(ne - 1)] for i = 1:model.dim.c]...)]
end
#
# function residual(model::ContactDynamicsModel, z, θ, κ)
# 	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
# 	q2, γ1, b1, η1, s1, s2 = unpack_z(model, z)
#
# 	ϕ = ϕ_func(model, q2)
# 	k = kinematics(model, q2)
#
# 	λ1 = [b1; γ1]
# 	vT = (J_fast(model, q2) * (q2 - q1) / h[1])[1:2]
#
# 	ne = dim(model.env)
#
# 	[model.dyn.d(h, q0, q1, u1, w1, λ1, q2);
# 	 vT - η1[2:3];
# 	 s1 - ϕ;
# 	 s2 - μ[1] * γ1;
# 	 γ1 .* s1 .- κ;
# 	 second_order_cone_product(η1, [s2; b1]) - [κ; 0.0; 0.0]]
# end

model = particle_nonlinear
num_var(model)
num_data(model)
idx_ineq = inequality_indices(model)
idx_soc = soc_indices(model)

q0 = rand(model.dim.q)
γ0 = rand(model.dim.c)
b0 = rand(model.dim.b)
η0 = rand(model.dim.b + model.dim.c)

z0 = zeros(num_var(model))
θ0 = zeros(num_data(model))

z_initialize!(z0, model, q0)

ne = dim(model.env)
bb = [z0[idx] for idx in idx_soc]

contact_forces(model, 10.0, [1.0; 2.0], q0, [0.0; 0.0])
velocity_stack(model, q0, rand(model.dim.q), [0.0; 0.0], 0.1)

dir = joinpath(@__DIR__, "particle_nonlinear_cone")
# model = deepcopy(particle_nonlinear)
#
# path_base = joinpath(dir, "dynamics/base.jld2")
# path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
# path_res = joinpath(dir, "flat/residual.jld2")
# path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
# path_linearized = joinpath(dir, "flat/linearized.jld2")

model = deepcopy(particle_nonlinear_quadratic)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "quadratic/residual.jld2")
path_jac = joinpath(dir, "quadratic/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "quadratic/linearized.jld2")

expr_base = generate_base_expressions(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

model.spa.rz_sp = rz_sp
model.spa.rθ_sp = rθ_sp

expr_linearized = generate_linearized_expressions(model)
save_expressions(expr_linearized, path_linearized, overwrite=true)
instantiate_linearized!(model, path_linearized)

# time
h = 0.01
T = 1000

## DROP
# initial conditions
q1 = @SVector [0.05, 1.0, 2.0]
q0 = @SVector [0.0, 1.0, 2.0]
ϕ_func(model, q1)
# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-5,
		diff_sol = true),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test status

plot(hcat(sim.traj.q[1:3:T]...)', label = ["x" "y" "z"],
	legend = :topright)
@show ϕ_func(model, sim.traj.q[end])
@show model.env.surf(sim.traj.q[end])
@show sim.traj.q[end]

include(joinpath(pwd(), "src/dynamics/particle_nonlinear_cone/visuals.jl"))

vis = Visualizer()
render(vis)
visualize!(vis, model, sim.traj.q,
	Δt = h, r = 0.1)
plot_surface!(vis, model.env,
	col = (1.0, 0.0, 0.0), α = 0.25,
	xlims = [-5.0, 5.0],
	ylims = [-5.0, 5.0])
