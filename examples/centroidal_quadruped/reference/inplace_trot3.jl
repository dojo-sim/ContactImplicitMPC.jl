# ## model
include("trajopt_model_v2.jl")

# ## horizon
h = 0.063636363
T = 12
Tm = 3 # mid point for a swing / stance change 
0.7 / 7
0.7 / 11
# ## centroidal_quadruped
s = get_simulation("centroidal_quadruped", "flat_3D_lc", "flat")
model = s.model
env = s.env

nx = 2 * model.nq
nc = 4 #model.nc
nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
nθ = 53

# ## model
d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dyn1(model, env, [h], y, x, u, w), nx + nθ + nx, nx, nu)
dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(model, env, [h], y, x, u, w), nx + nθ + nx, nx + nθ + nx, nu)

dyn = [d1, [dt for t = 2:T-1]...]

# ## initial conditions
body_height = 0.25
foot_x = 0.17
foot_y = 0.17
foot_height = 0.08

q1 = zeros(model.nq)

# body position
q1[1:3] = [0.0; 0.0; body_height]
q1[4:6] = [0.0; 0.0; 0.0]

# foot1
q1[7:9] = [foot_x; foot_y; 0]

# foot2
q1[10:12] = [foot_x; -foot_y; 0]

# foot3
q1[13:15] = [-foot_x; foot_y; 0]

# foot4
q1[16:18] = [-foot_x; -foot_y; 0]

# qM1 FL and RR swing 
# qM2 FR and RL swing  
qM1 = deepcopy(q1)
qM2 = deepcopy(q1)
qM1[6 + 3] += foot_height # front left
qM1[15 + 3] += foot_height # back right
qM2[9 + 3] += foot_height # front right
qM2[12 + 3] += foot_height # back left

# Terminal Q 
qT = copy(q1)

vis = Visualizer() 
render(vis)
visualize!(vis, model, [qM2], Δt=h);

# q_ref = [q1, qM1, qT,  
#          q1, qM2, qM2, qT]

q_ref = [q1, linear_interpolation(q1, qM1, Tm)..., linear_interpolation(qM1, qT, Tm)..., 
         linear_interpolation(q1, qM2, Tm)..., linear_interpolation(qM2, qT, Tm)...]

for τ = 2:7
    if τ < 4.5 
        t = τ 
    else
        t = 8 - τ + 1
    end
    @show t
    @show q_ref[τ][9] = θ[1] + θ[2] * t + θ[3] * t^2 + θ[4] * t^3 + θ[5] * t^4 + θ[6] * t^5
    @show q_ref[τ][18] = θ[1] + θ[2] * t + θ[3] * t^2 + θ[4] * t^3 + θ[5] * t^4 + θ[6] * t^5
end

for τ = 8:12
    t = τ - 6
    if t >= 4.5 
        t = 8 - t + 1 
    end
    q_ref[τ][12] = θ[1] + θ[2] * t + θ[3] * t^2 + θ[4] * t^3
    q_ref[τ][15] = θ[1] + θ[2] * t + θ[3] * t^2 + θ[4] * t^3
end

# a + b * x + c * x^2 + d * x^3 = z
# b + c * 2x + d 3 x^2 = \partial z

# x = 1, z = 0, dz = 0
# x = 4, z = 0.8, dz = 0

# A = [
# 	1.0 2.0 2.0^2 2.0^3;
# 	0.0 1.0 (2.0 * 2.0) (3.0 * 2.0^2);
# 	1.0 5.0 (5.0^2) (5.0^3)
# 	0.0 1.0 (2.0 * 5.0) (3.0 * 5.0^2);
# ]
# b = [
# 	0.0;
# 	0.0;
# 	0.08; 
# 	0.0;
# ]

A = [
	1.0 2.0 2.0^2 2.0^3 2.0^4 2.0^5;
    1.0 3.0 3.0^2 3.0^3 3.0^4 3.0^5;
	1.0 4.0 4.0^2 4.0^3 4.0^4 4.0^5;
	1.0 5.0 5.0^2 5.0^3 5.0^4 5.0^5;
    1.0 6.0 6.0^2 6.0^3 6.0^4 6.0^5;
	1.0 7.0 7.0^2 7.0^3 7.0^4 7.0^5;
]
b = [
	0.0;
    0.04
	0.08;
	0.08;
    0.04
	0.0;
]

rank(A)
θ = A \ b

# x = 1
# z = θ[1] + θ[2] * x + θ[3] * x^2 + θ[4] * x^3

# x = 2
# z = θ[1] + θ[2] * x + θ[3] * x^2 + θ[4] * x^3

# x = 3
# z = θ[1] + θ[2] * x + θ[3] * x^2 + θ[4] * x^3

# x = 3.5
# z = θ[1] + θ[2] * x + θ[3] * x^2 + θ[4] * x^3

# x = 4.5
# z = θ[1] + θ[2] * x + θ[3] * x^2 + θ[4] * x^3


x_ref = [[q_ref[t]; q_ref[t+1]] for t = 1:T]
x1 = x_ref[1]
xM = x_ref[Tm]
xT = x_ref[T]

visualize!(vis, model, q_ref, Δt=h);

plot(hcat(q_ref...)', labels="")

# ## objective
obj = DTO.Cost{Float64}[]
for t = 1:T
    if t == T
        function objT(x, u, w)
            J = 0.0
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-3 * dot(v, v)
            J += 100 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t])
            return J
        end
        push!(obj, DTO.Cost(objT, nx + nθ + nx, 0))
    elseif t == 1
        function obj1(x, u, w)
            J = 0.0
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-3 * dot(v, v)
            J += 100 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t])
            J += 0.5 * transpose(u[1:model.nu]) * Diagonal(1.0e-3 * ones(model.nu)) * u[1:model.nu]
            J += 1000.0 * u[end] # slack
            return J
        end
        push!(obj, DTO.Cost(obj1, nx, nu))
    else
        function objt(x, u, w)
            J = 0.0
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-3 * dot(v, v)
            u_previous = x[nx .+ (1:53)]
            u_control = u
            w = (u_control - u_previous) ./ h
            J += 0.5 * 1.0e-3 * dot(w, w)
            J += 100 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t])
            J += 0.5 * transpose(u[1:model.nu]) * Diagonal(1.0e-3 * ones(model.nu)) * u[1:model.nu]
            J += 1000.0 * u[end] # slack
            return J
        end
        push!(obj, DTO.Cost(objt, nx + nθ + nx, nu))
    end
end

# ## constraints
# initial condition
xl1 = [q1; q1]
xu1 = [q1; q1]
xlt = [-Inf * ones(nx); -Inf * ones(nθ); -Inf * ones(nx)]
xut = [Inf * ones(nx); Inf * ones(nθ); Inf * ones(nx)]

# final condition
xlT = [-Inf * ones(nq); qT; -Inf * ones(nθ); -Inf * ones(nx)]
xuT = [Inf * ones(nq); qT; Inf * ones(nθ); Inf * ones(nx)]

ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

# bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
# bndt = DTO.Bound(nx + nθ + nx, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
# bndT = DTO.Bound(nx + nθ + nx, 0, xl=xlT, xu=xuT)
# bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];

bnds = DTO.Bound{Float64}[]
push!(bnds, DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu))
for t = 2:T-1 
    push!(bnds, DTO.Bound(nx + nθ + nx, nu, 
        xl=[-Inf * ones(nq); -Inf * ones(nq); -Inf * ones(nθ); -Inf * ones(nx)], 
        xu=[Inf * ones(nq); Inf * ones(nq); Inf * ones(nθ); Inf * ones(nx)], 
        ul=ul, uu=uu))
end
push!(bnds, DTO.Bound(nx + nθ + nx, 0, xl=xlT, xu=xuT))


cons = DTO.Constraint{Float64}[]
for t = 1:T
    if t == 1
        function constraints_1(x, u, w)
            [
            # equality (16)
            contact_constraints_equality(model, env, h, x, u, w);
            # inequality (28)
            contact_constraints_inequality_1(model, env, h, x, u, w);

            # body/feet constraints
            x[3] - x_ref[t][3]; # body height
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            # x[9:3:18] - x_ref[t][9:3:18];
            ]
        end
        push!(cons, DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(16 .+ (1:28))))
    elseif t == T
        function constraints_T(x, u, w)
            [
            # inequality (8)
            contact_constraints_inequality_T(model, env, h, x, u, w);

            # body/feet constraints
            x[3] - x_ref[t][3]; # body height
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            # x[9:3:18] - x_ref[t][9:3:18];
            ]
        end
        push!(cons, DTO.Constraint(constraints_T, nx + nθ + nx, nu, idx_ineq=collect(0 .+ (1:8))))
    else
        function constraints_t(x, u, w)
            [
            # equality (16)
            contact_constraints_equality(model, env, h, x, u, w);
            # inequality (32)
            contact_constraints_inequality_t(model, env, h, x, u, w);

            # body/feet constraints
            x[3] - x_ref[t][3]; # body height
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            # x[9:3:18] - x_ref[t][9:3:18];
            ]
        end
        push!(cons, DTO.Constraint(constraints_t, nx + nθ + nx, nu, idx_ineq=collect(16 .+ (1:32))) )
    end
end

# ## problem
tolerance = 1.0e-3
p = DTO.solver(dyn, obj, cons, bnds,
    options=DTO.Options(
        max_iter=2000,
        tol=tolerance,
        constr_viol_tol=tolerance,
        ))

# ## initialize
x_interpolation = [x_ref[1], [[x_ref[t]; zeros(nθ); zeros(nx)] for t = 2:T]...]
u_guess = [1.0e-4 * rand(nu) for t = 1:T-1] # may need to run more than once to get good trajectory
DTO.initialize_states!(p, x_interpolation)
DTO.initialize_controls!(p, u_guess)

# ## solve
@time DTO.solve!(p)

# ## solution
x_sol, u_sol = DTO.get_trajectory(p)
@show x_sol[1]
@show x_sol[T]
maximum([u[end] for u in u_sol[1:end-1]])

# ## visualize
vis = Visualizer()
render(vis)
visualize!(vis, model, [x_sol[1][1:nq], [x[nq .+ (1:nq)] for x in x_sol]...], Δt=h);

q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
u_opt = [u[1:model.nu] for u in u_sol]
λ_opt = [u[model.nu .+ (1:4)] for u in u_sol]
b_opt = [u[model.nu + 4 .+ (1:16)] for u in u_sol]

q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
u_opt = [u[1:model.nu] for u in u_sol]
γ_opt = [u[model.nu .+ (1:4)] for u in u_sol]
b_opt = [u[model.nu + 4 .+ (1:16)] for u in u_sol]
ψ_opt = [u[model.nu + 4 + 16 .+ (1:4)] for u in u_sol]
η_opt = [u[model.nu + 4 + 16 + 4 .+ (1:16)] for u in u_sol]

qm = copy(q_opt)
um = copy(u_opt)
γm = copy(γ_opt)
bm = copy(b_opt)
ψm = copy(ψ_opt)
ηm = copy(η_opt)
μm = model.μ_world

hm = h
timesteps = range(0.0, stop=(h * (length(qm) - 2)), length=(length(qm) - 2))
plot(hcat(qm...)', labels="")

plot(timesteps, hcat(qm[2:end-1]...)', labels="")
plot(timesteps, hcat(um...)', labels="")
plot(timesteps, hcat(γm...)', labels="")
plot(timesteps, hcat(bm...)', labels="")
plot(timesteps, hcat(ψm...)', labels="")
plot(timesteps, hcat(ηm...)', labels="")

using JLD2
@save joinpath(@__DIR__, "inplace_trot_v7.jld2") qm um γm bm ψm ηm μm hm
@load joinpath(@__DIR__, "inplace_trot_v7.jld2") qm um γm bm ψm ηm μm hm
