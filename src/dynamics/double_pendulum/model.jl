"""
    Double pendulum
"""

struct DoublePendulum{T} <: ContactModel
    dim::Dimensions

    m1::T    # mass link 1
    J1::T    # inertia link 1
    l1::T    # length link 1
    lc1::T   # length to COM link 1

    m2::T    # mass link 2
    J2::T    # inertia link 2
    l2::T    # length link 2
    lc2::T   # length to COM link 2

    g::T     # gravity

    b1::T    # joint friction
    b2::T
end

lagrangian(model::DoublePendulum, q, q̇) = 0.0

function kinematics(model::DoublePendulum, x)
    @SVector [model.l1 * sin(x[1]) + model.l2 * sin(x[1] + x[2]),
              -1.0 * model.l1 * cos(x[1]) - model.l2 * cos(x[1] + x[2])]
end

function kinematics_elbow(model::DoublePendulum, x)
    @SVector [model.l1 * sin(x[1]),
              -1.0 * model.l1 * cos(x[1])]
end

function M_func(model::DoublePendulum, x)
    a = (model.J1 + model.J2 + model.m2 * model.l1 * model.l1
         + 2.0 * model.m2 * model.l1 * model.lc2 * cos(x[2]))

    b = model.J2 + model.m2 * model.l1 * model.lc2 * cos(x[2])

    c = model.J2

    @SMatrix [a b;
              b c]
end

function τ_func(model::DoublePendulum, x)
    a = (-1.0 * model.m1 * model.g * model.lc1 * sin(x[1])
         - model.m2 * model.g * (model.l1 * sin(x[1])
         + model.lc2 * sin(x[1] + x[2])))

    b = -1.0 * model.m2 * model.g * model.lc2 * sin(x[1] + x[2])

    @SVector [a,
              b]
end

function c_func(model::DoublePendulum, q, q̇)
    a = -2.0 * model.m2 * model.l1 * model.lc2 * sin(q[2]) * q̇[2]
    b = -1.0 * model.m2 * model.l1 * model.lc2 * sin(q[2]) * q̇[2]
    c = model.m2 * model.l1 * model.lc2 * sin(q[2]) * q̇[1]
    d = 0.0

    @SMatrix [a b;
              c d]
end

function B_func(model::DoublePendulum, x)
    @SMatrix [0.0;
              1.0]
end

function C_func(model::DoublePendulum, q, q̇)
    c_func(model, q, q̇) * q̇ - τ_func(model, q)
end

function ϕ_func(model, q)
    # SVector{model.dim.c}([q[1], 0.5 * π - q[2], q[2] + 0.5 * π])
    SVector{model.dim.c}([0.5 * π - q[2], q[2] + 0.5 * π])
end

function P_func(model, q)
    ϕ(z) = ϕ_func(model, z)
    ForwardDiff.jacobian(ϕ, q)
end

nq = 2
nu = 1
nw = 0
nc = 1
nb = 0

model = DoublePendulum(Dimensions(nq, nu, nw, nc, nquat),
    1.0, 0.33, 1.0, 0.5, 1.0, 0.33, 1.0, 0.5, 9.81, 0.1, 0.1)
env = flat_2D_lc

s = Simulation(model, env)

function lagrangian_derivatives(model::DoublePendulum, q, v)
	D1L = -1.0 * C_func(model, q, v)
    D2L = M_func(model, q) * v
	return D1L, D2L
end

function dynamics(model::DoublePendulum, h, q0, q1, u1, λ1, q2)
	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	# return 0.0
	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ B_func(model, qm2) * u1
        + transpose(P_func(model, q2)) * λ1
        - h[1] * 0.5 .* vm2)
end

function residual(model::DoublePendulum, env::Environment{<:World,LinearizedCone}, z, θ, κ)
	# nc = model.dim.c
	# nb = nc * friction_dim(env)
	# nf = Int(nb / nc)
	# np = dim(env)
    nq = model.dim.q
    nu = model.dim.u
    nc = model.dim.c

    q0 = θ[1:nq]
    q1 = θ[nq .+ (1:nq)]
    u1 = θ[2nq .+ (1:nu)]
    h = θ[2nq + nu .+ (1:1)]

    q2 = z[1:nq]
    λ1 = z[nq .+ (1:nc)]
    s1 = z[nq + nc .+ (1:nc)]

	# q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	# q2, γ1, b1, ψ1, s1, η1, s2 = unpack_z(model, env, z)

	# ϕ = ϕ_func(model, env, q2)
    #
	# k = kinematics(model, q2)
	# λ1 = contact_forces(model, env, γ1, b1, q2, k)
	# Λ1 = transpose(J_func(model, env, q2)) * λ1 #@@@@ maybe need to use J_fast
	# vT_stack = velocity_stack(model, env, q1, q2, k, h)
	# ψ_stack = transpose(E_func(model, env)) * ψ1
    # return dynamics(model, h, q0, q1, u1, zeros(model.dim.c), q2)

    [
     dynamics(model, h, q0, q1, u1, λ1, q2);
     s1 .- ϕ_func(model, q2);
     λ1 .* s1 .- κ;
    ]
	# [
    #  dynamics(model, h, q0, q1, u1, q2);
	#  # s1 - ϕ;
	#  # vT_stack + ψ_stack - η1;
	#  # s2 .- (μ[1] * γ1 .- E_func(model, env) * b1);
	#  # γ1 .* s1 .- κ;
	#  # b1 .* η1 .- κ;
	#  # ψ1 .* s2 .- κ
    #  ]
end

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nquat = 0
nz = nq + nc + nc
nθ = nq + nq + nu + 1

# Declare variables
@variables z[1:nz]
@variables θ[1:nθ]
@variables κ[1:1]

# Residual
r = residual(model, env, z, θ, κ)
r = Symbolics.simplify.(r)
rz = Symbolics.jacobian(r, z, simplify = true)
rθ = Symbolics.jacobian(r, θ, simplify = true)

# Build function
r_func = eval(build_function(r, z, θ, κ)[2])
rz_func = eval(build_function(rz, z, θ)[2])
rθ_func = eval(build_function(rθ, z, θ)[2])

q0 = [0.5 * π; 0.0]
q1 = [0.5 * π; 0.0]
u1 = [0.0]
h = 0.1

z0 = copy([q1; ones(nc + nc)])
θ0 = [q0; q1; u1; h]

# options
opts = InteriorPointOptions(
   κ_init = 0.1,
   κ_tol = 1.0e-4,
   r_tol = 1.0e-8,
   diff_sol = true)

idx_ineq = collect([nq .+ (1:(nc + nc))]...)

# solver
ip = interior_point(z0, θ0,
   r! = r_func, rz! = rz_func,
   rz = similar(rz, Float64),
   rθ! = rθ_func,
   rθ = similar(rθ, Float64),
   idx_ineq = idx_ineq,
   opts = opts)

s.res.r! = r_func
s.res.rz! = rz_func
s.res.rθ! = rθ_func
s.rz = ip.rz
s.rθ = ip.rθ

# # r_func(ones(nz), z0, θ0, [1.0])
# # rz_func(ones(nz, nz), z0, θ0)
# # rθ_func(ones(nz, nθ), z0, θ0)
#
#simulate
# T = 50
# q_hist = [q0, q1]
#
# for t = 1:T
#     ip.z .= copy([q_hist[end]; 0.1 * ones(nc + nc)])
#     ip.θ .= copy([q_hist[end-1]; q_hist[end]; u1; h])
#     status = interior_point_solve!(ip)
#
#     if status
#         push!(q_hist, ip.z[1:nq])
#     else
#         println("dynamics failure t = $t")
#         println("res norm = $(norm(ip.r, Inf))")
#         break
#     end
# end
#
# # visualize!(vis, model, q_hist, Δt = h)
#
#
# vis = Visualizer()
# render(vis)
# default_background!(vis)
# settransform!(vis["/Cameras/default"],
#         compose(Translation(0.0, -95.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
# setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 30)
#
# visualize!(vis, model, q_hist, Δt = h)
