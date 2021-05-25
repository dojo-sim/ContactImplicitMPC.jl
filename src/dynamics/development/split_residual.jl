model = deepcopy(hopper_2D_sinusoidal)
dir = joinpath(pwd(), "src/dynamics/hopper_2D")

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "sinusoidal/residual.jld2")
path_jac = joinpath(dir, "sinusoidal/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

expr_base = generate_base_expressions_analytical(model)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model, derivs = true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn, derivs = true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)
model.spa.rz_sp = rz_sp
model.spa.rθ_sp = rθ_sp

expr_linearized = generate_linearized_expressions(model)
save_expressions(expr_linearized, path_linearized, overwrite=true)
instantiate_linearized!(model, path_linearized)

cf_expr = generate_contact_expressions(model)
cf_methods = ContactMethods()
instantiate_contact_methods!(cf_methods, cf_expr)

_z = rand(num_var(model))
_θ = rand(num_data(model))
_q0, _q1, _u1, _w1, _μ, _h = unpack_θ(model, _θ)
_q2, _γ1, _b1, _ψ1, _η1, _s1, _s2 = unpack_z(model, _z)
_λ1 = cf_methods.cf(_γ1, _b1, _q2)
_vT = cf_methods.vs(_q1, _q2, _h)
model.res.rz!(model.spa.rz_sp, _z, _θ)
model.spa.rz_sp
rank(model.spa.rz_sp)

model.res.rθ!(model.spa.rθ_sp, _z, _θ)
model.spa.rθ_sp
rank(model.spa.rθ_sp)

nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c
nb = model.dim.b
ncf = nc * dim(model.env)

# Declare variables
# @variables q0[1:nq]
# @variables q1[1:nq]
# @variables u1[1:nu]
# @variables w1[1:nw]
# @variables γ1[1:nc]
# @variables λ1[1:ncf]
# @variables b1[1:nb]
# @variables q2[1:nq]
# @variables vs[1:ncf]
# @variables ψ1[1:nc]
# @variables η1[1:nb]
# @variables h[1:1]

# Dynamics
# d = dynamics(model, h, q0, q1, u1, w1, γ1, b1, q2)
# d = dynamics(model, h, q0, q1, u1, w1, λ1, q2)
# d = Symbolics.simplify.(d)
# dy  = Symbolics.jacobian(d, [q0; q1; u1; w1; λ1; q2], simplify=true)
# dq2 = Symbolics.jacobian(d, q2, simplify = true)
# dλ1 = Symbolics.jacobian(d, λ1, simplify = true)
# dθ = Symbolics.jacobian(d, [q0; q1; u1; w1; h], simplify = true)
#
# dq2_func  = eval(build_function(dq2, h, q0, q1, u1, w1, λ1, q2)[1])
# dλ1_func  = eval(build_function(dλ1, h, q0, q1, u1, w1, λ1, q2)[1])
# dθ_func   = eval(build_function(dθ, h, q0, q1, u1, w1, λ1, q2)[1])
#
# cf = contact_forces(model, γ1, b1, q2)
# dcf  = Symbolics.jacobian(cf, [q2; γ1; b1], simplify=true)
# dcf_func = eval(build_function(dcf, γ1, b1, q2)[1])

# _vs = velocity_stack(model, q1, q2, h)
# vsq2 = Symbolics.jacobian(_vs, q2, simplify = true)
# vsq1h = Symbolics.jacobian(_vs, [q1; h], simplify = true)
# vsq2_func = eval(build_function(vsq2, q1, q2, h)[1])
# vsq1h_func = eval(build_function(vsq1h, q1, q2, h)[1])
#
# dq2_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2)
# dλ1_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2) * dcf_func(_q2, _γ1, _b1)
# dθ_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2)

rd = zeros(model.dim.q, num_var(model))
rd[:, 1:model.dim.q] = model.dyn.dq2(_h, _q0, _q1, _u1, _w1, _λ1, _q2)
rd[:, 1:(model.dim.q + model.dim.c + model.dim.b)] += model.dyn.dλ1(_h, _q0, _q1, _u1, _w1, _λ1, _q2) * cf_methods.dcf(_γ1, _b1,_q2)
@test norm(model.spa.rz_sp[1:model.dim.q, :] - rd) < 1.0e-8

# md = vs + transpose(E_func(model)) * ψ1 - η1
# mdvs = Symbolics.jacobian(md, vs, simplify = true)
# mdvs_func = eval(build_function(mdvs, vs, ψ1, η1)[1])
# mdψη = Symbolics.jacobian(md, [ψ1; η1], simplify = true)
# mdψη_func = eval(build_function(mdψη, vs, ψ1, η1)[1])

rmd = zeros(model.dim.b, num_var(model))
rmd[:, model.dim.q + model.dim.c + model.dim.b .+ (1:model.dim.c + model.dim.b)] = cf_methods.mdψη(_vT, _ψ1, _η1)# copy(mdψη_func(_vT, _ψ1, _η1))
rmd[:, 1:model.dim.q] += cf_methods.mdvs(_vT, _ψ1, _η1) * cf_methods.vsq2(_q1, _q2, _h)

@test norm(model.spa.rz_sp[model.dim.q .+ (1:model.dim.b), :] - rmd) < 1.0e-8

rθd = zeros(model.dim.q, num_data(model))
idx = collect([(1:2model.dim.q + model.dim.u + model.dim.w)..., num_data(model)])
rθd[:, idx] = model.dyn.dθ(_h, _q0, _q1, _u1, _w1, _λ1, _q2)

@test norm(model.spa.rθ_sp[1:model.dim.q, :] - rθd[:, :]) < 1.0e-8

rθmd = zeros(model.dim.b, num_data(model))
# vsq1h_func(_q1, _q2, _h)
idx = collect([(model.dim.q .+ (1:model.dim.q))..., (2model.dim.q + model.dim.u + model.dim.w + 1 .+ (1:1))...])
rθmd[:, idx] = cf_methods.vsq1h(_q1, _q2, _h)

@test norm(model.spa.rθ_sp[model.dim.q .+ (1:model.dim.b), :] - rθmd) < 1.0e-8

# function res_con(model::ContactModel, z, θ, κ)
# 	nc = model.dim.c
# 	nb = model.dim.b
# 	nf = Int(nb / nc)
# 	np = dim(model.env)
#
# 	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
# 	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)
#
# 	[s1 - ϕ_func(model, q2);
# 	 s2 .- (μ[1] * γ1 .- E_func(model) * b1);
# 	 γ1 .* s1 .- κ;
# 	 b1 .* η1 .- κ;
# 	 ψ1 .* s2 .- κ]
# end

nz = num_var(model)
nθ = num_data(model)

# # Declare variables
# @variables z[1:nz]
# @variables θ[1:nθ]
# @variables κ

# Residual
# rc = res_con(model, z, θ, κ)
# rc = Symbolics.simplify.(rc)
# rcz = Symbolics.jacobian(rc, z, simplify = true)
# rcθ = Symbolics.jacobian(rc, θ, simplify = true) # TODO: sparse version
#
# rcz_sp = similar(rcz, Float64)
# rcθ_sp = similar(rcθ, Float64)

# rc_func = eval(build_function(rc, z, θ, κ)[2])
# rcz_func = eval(build_function(rcz, z, θ)[2])
# rcθ_func = eval(build_function(rcθ, z, θ)[2])
#
# rcz_func(rcz_sp, _z, _θ)
# rcθ_func(rcθ_sp, _z, _θ)

rcz_sp = zeros(nz - model.dim.q - model.dim.b, nz)
rcθ_sp = zeros(nz - model.dim.q - model.dim.b, nθ)

cf_methods.rcz(rcz_sp, _z, _θ)
cf_methods.rcθ(rcθ_sp, _z, _θ)

@test norm(model.spa.rz_sp[(model.dim.q + model.dim.b + 1):end, :] - rcz_sp) < 1.0e-8
@test norm(model.spa.rθ_sp[(model.dim.q + model.dim.b + 1):end, :] - rcθ_sp) < 1.0e-8

function rz_split!(rz, z, θ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)
	λ1 = cf_methods.cf(γ1, b1, q2)
	vT = cf_methods.vs(q1, q2, h)


	# Dynamics
	rz[1:model.dim.q, 1:model.dim.q] = model.dyn.dq2(h, q0, q1, u1, w1, λ1, q2)
	rz[1:model.dim.q, 1:(model.dim.q + model.dim.c + model.dim.b)] += model.dyn.dλ1(h, q0, q1, u1, w1, λ1, q2) * cf_methods.dcf(γ1, b1, q2)

	# Maximum dissipation
	rz[model.dim.q .+ (1:model.dim.b), model.dim.q + model.dim.c + model.dim.b .+ (1:model.dim.c + model.dim.b)] = cf_methods.mdψη(vT, ψ1, η1)# copy(mdψη_func(_vT, _ψ1, _η1))
	rz[model.dim.q .+ (1:model.dim.b), 1:model.dim.q] += cf_methods.mdvs(vT, ψ1, η1) * cf_methods.vsq2(q1, q2, h)

	# Other constraints
	cf_methods.rcz(view(rz, (model.dim.q + model.dim.b + 1):num_var(model), :), z, θ)
end

function rθ_split!(rθ, z, θ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)
	λ1 = cf_methods.cf(γ1, b1, q2)
	vT = cf_methods.vs(q1, q2, h)

	# Dynamics
	idx = collect([(1:2model.dim.q + model.dim.u + model.dim.w)..., num_data(model)])
	rθ[1:model.dim.q, idx] = model.dyn.dθ(h, q0, q1, u1, w1, λ1, q2)

	# Maximum dissipation
	idx = collect([(model.dim.q .+ (1:model.dim.q))..., (2model.dim.q + model.dim.u + model.dim.w + 1 .+ (1:1))...])
	rθ[model.dim.q .+ (1:model.dim.b), idx] = cf_methods.vsq1h(q1, q2, h)

	# Other constraints
	cf_methods.rcθ(view(rθ, (model.dim.q + model.dim.b + 1):num_var(model), :), z, θ)
end

rz_test = similar(model.spa.rz_sp, Float64)
rθ_test = similar(model.spa.rθ_sp, Float64)

rz_split!(rz_test, _z, _θ)
rθ_split!(rθ_test, _z, _θ)

@test norm(model.spa.rz_sp - rz_test) < 1.0e-8
@test norm(model.spa.rθ_sp - rθ_test) < 1.0e-8
