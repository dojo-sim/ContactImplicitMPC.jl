model = deepcopy(quadruped)
dir = joinpath(pwd(), "src/dynamics/quadruped")
model = get_model("quadruped", surf="flat")
include(joinpath(pwd(), "src/simulator/environment/piecewise.jl"))
model.env = Environment{R2}(terrain_sym, d_terrain_sym)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "piecewise/residual.jld2")
path_jac = joinpath(dir, "piecewise/sparse_jacobians.jld2")
# path_linearized = joinpath(dir, "sinusoidal/linearized.jld2")

instantiate_base!(model, path_base)
instantiate_dynamics!(model, path_dyn)
# expr_base = generate_base_expressions_analytical(model)
# save_expressions(expr_base, path_base, overwrite=true)
# instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model, derivs = true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn, derivs = true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model, jacobians = false)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res, jacobians = false)
model.spa.rz_sp = rz_sp
model.spa.rθ_sp = rθ_sp

# expr_linearized = generate_linearized_expressions(model)
# save_expressions(expr_linearized, path_linearized, overwrite=true)
# instantiate_linearized!(model, path_linearized)

cf_expr = generate_contact_expressions(model)
cf_methods = ContactMethods()
instantiate_contact_methods!(cf_methods, cf_expr)

_z = rand(num_var(model))
_θ = rand(num_data(model))
_q0, _q1, _u1, _w1, _μ, _h = unpack_θ(model, _θ)
_q2, _γ1, _b1, _ψ1, _η1, _s1, _s2 = unpack_z(model, _z)
_λ1 = cf_methods.cf(_γ1, _b1, _q2)
_vT = cf_methods.vs(_q1, _q2, _h)

nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c
nb = model.dim.b
ncf = nc * dim(model.env)

rd = zeros(model.dim.q, num_var(model))
rd[:, 1:model.dim.q] = model.dyn.dq2(_h, _q0, _q1, _u1, _w1, _λ1, _q2)
rd[:, 1:(model.dim.q + model.dim.c + model.dim.b)] += model.dyn.dλ1(_h, _q0, _q1, _u1, _w1, _λ1, _q2) * cf_methods.dcfa(_γ1, _b1, _q2)

rmd = zeros(model.dim.b, num_var(model))
rmd[:, model.dim.q + model.dim.c + model.dim.b .+ (1:model.dim.c + model.dim.b)] = cf_methods.mdψη(_vT, _ψ1, _η1)# copy(mdψη_func(_vT, _ψ1, _η1))
rmd[:, 1:model.dim.q] += cf_methods.mdvs(_vT, _ψ1, _η1) * cf_methods.vsaq2(_q1, _q2, _h) # NOTE: breaks?

rθd = zeros(model.dim.q, num_data(model))
idx = collect([(1:2model.dim.q + model.dim.u + model.dim.w)..., num_data(model)])
rθd[:, idx] = model.dyn.dθ(_h, _q0, _q1, _u1, _w1, _λ1, _q2)

rθmd = zeros(model.dim.b, num_data(model))
idx = collect([(model.dim.q .+ (1:model.dim.q))..., (2model.dim.q + model.dim.u + model.dim.w + 1 .+ (1:1))...])
rθmd[:, idx] = cf_methods.vsaq1h(_q1, _q2, _h)

nz = num_var(model)
nθ = num_data(model)

# # Declare variables
@variables z[1:nz]
@variables θ[1:nθ]
@variables κ[1:1]

# Residual

rc = res_con(model, z, θ, κ)
rc = Symbolics.simplify.(rc)
rcz = Symbolics.jacobian(rc, z, simplify = true)
rcθ = Symbolics.jacobian(rc, θ, simplify = true) # TODO: sparse version
#
rcz_sp = similar(rcz, Float64)
rcθ_sp = similar(rcθ, Float64)

rc_func = eval(build_function(rc, z, θ, κ)[2])
rcz_func = eval(build_function(rcz, z, θ)[2])
rcθ_func = eval(build_function(rcθ, z, θ)[2])

rcz_func(rcz_sp, _z, _θ)
rcθ_func(rcθ_sp, _z, _θ)

rcz_sp = zeros(nz - model.dim.q - model.dim.b, nz)
rcθ_sp = zeros(nz - model.dim.q - model.dim.b, nθ)

cf_methods.rcz(rcz_sp, _z, _θ)
cf_methods.rcθ(rcθ_sp, _z, _θ)

# @test norm(model.spa.rz_sp[(model.dim.q + model.dim.b + 1):end, :] - rcz_sp) < 1.0e-8
# @test norm(model.spa.rθ_sp[(model.dim.q + model.dim.b + 1):end, :] - rcθ_sp) < 1.0e-8

function rz_split!(rz, z, θ)
	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)
	λ1 = cf_methods.cf(γ1, b1, q2)
	vT = cf_methods.vs(q1, q2, h)


	# Dynamics
	rz[1:model.dim.q, 1:model.dim.q] = model.dyn.dq2(h, q0, q1, u1, w1, λ1, q2)
	rz[1:model.dim.q, 1:(model.dim.q + model.dim.c + model.dim.b)] += model.dyn.dλ1(h, q0, q1, u1, w1, λ1, q2) * cf_methods.dcfa(γ1, b1, q2)

	# Maximum dissipation
	rz[model.dim.q .+ (1:model.dim.b), model.dim.q + model.dim.c + model.dim.b .+ (1:model.dim.c + model.dim.b)] = cf_methods.mdψη(vT, ψ1, η1)# copy(mdψη_func(_vT, _ψ1, _η1))
	rz[model.dim.q .+ (1:model.dim.b), 1:model.dim.q] += cf_methods.mdvs(vT, ψ1, η1) * cf_methods.vsaq2(q1, q2, h)

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
	rθ[model.dim.q .+ (1:model.dim.b), idx] = cf_methods.vsaq1h(q1, q2, h)

	# Other constraints
	cf_methods.rcθ(view(rθ, (model.dim.q + model.dim.b + 1):num_var(model), :), z, θ)
end

rz_test = similar(model.spa.rz_sp, Float64)
rz_test .= 0.0
rθ_test = similar(model.spa.rθ_sp, Float64)
rθ_test .= 0.0

rz_split!(rz_test, _z, _θ)
rθ_split!(rθ_test, _z, _θ)

# @test norm(model.spa.rz_sp[10:end,:] - rz_test[10:end,:]) < 1.0e-8
# @test norm(model.spa.rθ_sp - rθ_test) < 1.0e-8
