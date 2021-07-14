using BenchmarkTools, InteractiveUtils

include(joinpath(module_dir(), "src", "dynamics", "quadruped", "model.jl"))

s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env

q = rand(model.dim.q)
M = zeros(model.dim.q, model.dim.q)

@benchmark $M .= M_func($model, $q)


function test(x)
	transpose(x) * x
end

n = 5
@variables x[1:n]

t = test(x)
dt = Symbolics.gradient(t, x, simplify=true)
ddt = Symbolics.hessian(t, x, simplify=true)

dt_f = eval(Symbolics.build_function(dt, x)[1])
ddt_f = eval(Symbolics.build_function(ddt, x)[1])

x0 = rand(n)

@benchmark dt_f($x0)
@benchmark ddt_f($x0)


r = zeros(num_var(model, env))
z = zeros(num_var(model, env))

θ = zeros(num_data(model))
κ = 1.0
@benchmark s.res.r!($r, $z, $θ, $κ)

@variables q[1:model.dim.q], q̇[1:model.dim.q], v[1:model.dim.q]
l = lagrangian(model, q, q̇)

dl = Symbolics.gradient(l, [q; q̇])
dlq̇ = dl[model.dim.q .+ (1:model.dim.q)]
dlq̇v = transpose(dlq̇) * v
ddlq̇v = Symbolics.gradient(dlq̇v, q̇)
Mv = eval(Symbolics.build_function(ddlq̇v, q, v)[1])

ddl = Symbolics.hessian(l, [q; q̇])
m = ddl[model.dim.q .+ (1:model.dim.q), model.dim.q .+ (1:model.dim.q)]

M = eval(Symbolics.build_function(m, q)[1])

q0 = rand(model.dim.q)
q̇0 = rand(model.dim.q)
v0 = rand(model.dim.q)

@benchmark M($q0) * $v0











@benchmark Mv($q0, $v0)














norm(Mv(q0, v0) - mf(q0) * v0)
