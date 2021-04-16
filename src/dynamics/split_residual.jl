model = get_model("quadruped", surf="sinusoidal")

_z = rand(num_var(model))
_θ = rand(num_data(model))
_q0, _q1, _u1, _w1, _μ, _h = unpack_θ(model, _θ)
_q2, _γ1, _b1, _ψ1, _η1, _s1, _s2 = unpack_z(model, _z)
_λ1 = contact_forces(model, _γ1, _b1, _q2)
_vT = velocity_stack(model, _q1, _q2, _h)
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
@variables q0[1:nq]
@variables q1[1:nq]
@variables u1[1:nu]
@variables w1[1:nw]
@variables γ1[1:nc]
@variables λ1[1:ncf]
@variables b1[1:nb]
@variables q2[1:nq]
@variables vs[1:ncf]
@variables ψ1[1:nc]
@variables η1[1:nb]
@variables h[1:1]

# Dynamics
# d = dynamics(model, h, q0, q1, u1, w1, γ1, b1, q2)
d = dynamics(model, h, q0, q1, u1, w1, λ1, q2)
d = Symbolics.simplify.(d)
# dy  = Symbolics.jacobian(d, [q0; q1; u1; w1; λ1; q2], simplify=true)
dq2 = Symbolics.jacobian(d, q2, simplify = true)
dλ1 = Symbolics.jacobian(d, λ1, simplify = true)
dθ = Symbolics.jacobian(d, [q0; q1; u1; w1; h], simplify = true)

dq2_func  = eval(build_function(dq2, h, q0, q1, u1, w1, λ1, q2)[1])
dλ1_func  = eval(build_function(dλ1, h, q0, q1, u1, w1, λ1, q2)[1])
dθ_func   = eval(build_function(dθ, h, q0, q1, u1, w1, λ1, q2)[1])

cf = contact_forces(model, γ1, b1, q2)
dcf  = Symbolics.jacobian(cf, [q2; γ1; b1], simplify=true)
dcf_func = eval(build_function(dcf, q2, γ1, b1)[1])

vt = velocity_stack(model, q1, q2, h)
vtq2 = Symbolics.jacobian(vt, q2, simplify = true)
vtq1h = Symbolics.jacobian(vt, [q1; h], simplify = true)
vtq2_func = eval(build_function(vtq2, q1, q2, h)[1])
vtq1h_func = eval(build_function(vtq1h, q1, q2, h)[1])

dq2_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2)
dλ1_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2) * dcf_func(_q2, _γ1, _b1)
dθ_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2)

rd = zeros(model.dim.q, num_var(model))
rd[:, 1:model.dim.q] = copy(dq2_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2))
rd[:, 1:(model.dim.q + model.dim.c + model.dim.b)] += copy(dλ1_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2) * dcf_func(_q2, _γ1, _b1))
norm(model.spa.rz_sp[1:model.dim.q, :] - rd)

function max_diss(vs, ψ, η)
	ψ_stack = transpose(E_func(model)) * ψ1
	vs + ψ_stack - η
end

md = max_diss(vs, ψ1, η1)
mdvs = Symbolics.jacobian(md, vs, simplify = true)
mdvs_func = eval(build_function(mdvs, vs, ψ1, η1)[1])
mdψη = Symbolics.jacobian(md, [ψ1; η1], simplify = true)
mdψη_func = eval(build_function(mdψη, vs, ψ1, η1)[1])

rmd = zeros(model.dim.b, num_var(model))
rmd[:, model.dim.q + model.dim.c + model.dim.b .+ (1:model.dim.c + model.dim.b)] = copy(mdψη_func(_vT, _ψ1, _η1))
rmd[:, 1:model.dim.q] += copy(mdvs_func(_vT, _ψ1, _η1) * vtq2_func(_q1, _q2, _h))

norm(model.spa.rz_sp[model.dim.q .+ (1:model.dim.b), :] - rmd)

rθd = zeros(model.dim.q, num_data(model))
idx = collect([(1:2model.dim.q + model.dim.u + model.dim.w)..., num_data(model)])
rθd[:, idx] = copy(dθ_func(_h, _q0, _q1, _u1, _w1, _λ1, _q2))

norm(model.spa.rθ_sp[1:model.dim.q, :] - rθd[:, :])

rθmd = zeros(model.dim.b, num_data(model))
vtq1h_func(_q1, _q2, _h)
idx = collect([(model.dim.q .+ (1:model.dim.q))..., (2model.dim.q + model.dim.u + model.dim.w + 1 .+ (1:1))...])
rθmd[:, idx] = copy(vtq1h_func(_q1, _q2, _h))

norm(model.spa.rθ_sp[model.dim.q .+ (1:model.dim.b), :] - rθmd)

function res_con(model::ContactDynamicsModel, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	np = dim(model.env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)

	[s1 - ϕ_func(model, q2);
	 s2 .- (μ[1] * γ1 .- E_func(model) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

nz = num_var(model)
nθ = num_data(model)

# Declare variables
@variables z[1:nz]
@variables θ[1:nθ]
@variables κ

# Residual
rc = res_con(model, z, θ, κ)
rc = Symbolics.simplify.(rc)
rcz = Symbolics.jacobian(rc, z, simplify = true)
rcθ = Symbolics.jacobian(rc, θ, simplify = true) # TODO: sparse version

rcz_sp = similar(rcz, T)
rcθ_sp = similar(rcθ, T)

# rc_func = eval(build_function(rc, z, θ, κ)[2])
rcz_func = eval(build_function(rcz, z, θ)[2])
rcθ_func = eval(build_function(rcθ, z, θ)[2])

rcz_func(rcz_sp, _z, _θ)
rcθ_func(rcθ_sp, _z, _θ)
