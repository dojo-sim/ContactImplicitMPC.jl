# ## model
include("trajopt_model_wall_4feet.jl")
include(joinpath(@__DIR__, "..", "..", "..", "src/dynamics/centroidal_quadruped_wall/visuals.jl"))
@load joinpath(@__DIR__, "wall_stand_FL_4.jld2") qm um γm bm ψm ηm μm hm

vis = Visualizer()
render(vis)

# ## horizon
T = 21
h = 0.05

# ## centroidal_quadruped
s = get_simulation("centroidal_quadruped_wall", "flat_3D_lc", "flat")
model = s.model
env = s.env
nx = 2 * model.nq
nc = model.nc
nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
nθ = 93

# ## model
d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(
	model, env, [h], y, x, u, w), nx + nθ, nx, nu)

dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(
	model, env, [h], y, x, u, w), nx + nθ , nx + nθ , nu)

dyn = [d1, [dt for t = 2:T-1]...]

# ## initial conditions
body_height = 0.3
foot_x = 0.17
foot_y = 0.17

function nominal_configuration(model::CentroidalQuadrupedWall)
    qm[end]
end

function final_configuration(model::CentroidalQuadrupedWall)
    x_shift = 0.0
    [
        qm[end][1:9];
        0.25;-foot_y; body_height;
       -foot_x + x_shift; foot_y; 0.0;
       -foot_x + x_shift;-foot_y; 0.0;
    ]
end


q1 = nominal_configuration(model)
qM = nominal_configuration(model)
qT = final_configuration(model)
q_ref = nominal_configuration(model)

visualize!(vis, model, [qT], Δt=h)

x1 = [q1; q1]
xM = [qM; qM]
xT = [qT; qT]
x_ref = [q_ref; q_ref]

q1[10]
rad = qT[10] - q1[10]
θ_rad = range(1.0 * π, stop = 0.5 * π, length=T)
foot_arc_x = rad * cos.(θ_rad) .+ q1[10] .+ rad
foot_arc_z = rad * sin.(θ_rad) .* body_height ./ rad# .+ q1[9]

qT[10]
plot(foot_arc_x, foot_arc_z, aspect_ratio=:equal)

# ## objective
obj = DTO.Cost{Float64}[]
for t = 1:T
    if t == 1
        push!(obj, DTO.Cost((x, u, w) -> begin
            J = 0.0
            J += 0.5 * transpose(x[1:6] - x_ref[1:6]) * Diagonal(10000.0 * ones(6)) * (x[1:6] - x_ref[1:6])

            J += 0.5 * 10000.0 * (x[10] - foot_arc_x[t])^2
            J += 0.5 * 10000.0 * (x[11] + foot_y)^2
            J += 0.5 * 10000.0 * (x[12] - foot_arc_z[t])^2

            J += 0.5 * transpose(u) * Diagonal([1.0 * ones(model.nu); zeros(nu - model.nu)]) * u
            J += 10000.0 * u[end] # slack
            J += 0.5 * transpose(u[model.nu + 8 .+ (1:40)]) * Diagonal(1.0e-3 * ones(40)) * u[model.nu + 8 .+ (1:40)]

            return J
        end,
        nx, nu))
    elseif t == T
        push!(obj, DTO.Cost((x, u, w) -> begin
            J = 0.0
            J += 0.5 * transpose(x[1:6] - xT[1:6]) * Diagonal(10000.0 * ones(6)) * (x[1:6] - xT[1:6])

            J += 0.5 * 10000.0 * (x[10] - foot_arc_x[t])^2
            J += 0.5 * 10000.0 * (x[11] + foot_y)^2
            J += 0.5 * 10000.0 * (x[12] - foot_arc_z[t])^2

            J -= 100.0 * x[36 + model.nu + 5]
            J -= 100.0 * x[36 + model.nu + 6]
            return J
        end, nx + nθ, 0))
    else
        push!(obj, DTO.Cost((x, u, w) -> begin
            J = 0.0

            u_prev = x[nx .+ (1:93)]
            w = (u - u_prev) ./ h
            J += 0.5 * 1.0 * dot(w[1:end-1], w[1:end-1])

            J += 0.5 * transpose(x[1:6] - x_ref[1:6]) * Diagonal(10000.0 * ones(6)) * (x[1:6] - x_ref[1:6])

            J += 0.5 * 10000.0 * (x[10] - foot_arc_x[t])^2
            J += 0.5 * 10000.0 * (x[11] + foot_y)^2
            J += 0.5 * 10000.0 * (x[12] - foot_arc_z[t])^2

            J += 0.5 * transpose(u) * Diagonal([1.0 * ones(model.nu); zeros(nu - model.nu)]) * u
            J += 10000.0 * u[end] # slack
            J += 0.5 * transpose(u[model.nu + 8 .+ (1:40)]) * Diagonal(1.0e-3 * ones(40)) * u[model.nu + 8 .+ (1:40)]

            return J
        end, nx + nθ, nu))
    end
end
# function obj1(x, u, w)
# 	J = 0.0
# 	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal([1000.0 * ones(6); 100.0 * ones(nx-6)]) * (x[1:nx] - x_ref)
# 	J += 0.5 * transpose(u) * Diagonal([1.0 * ones(model.nu); zeros(nu - model.nu)]) * u
#     J += 1000.0 * u[end] # slack
#     J += 0.5 * transpose(u[model.nu + 8 .+ (1:40)]) * Diagonal(1.0e-3 * ones(40)) * u[model.nu + 8 .+ (1:40)]

# 	return J
# end

# function objt(x, u, w)
# 	J = 0.0

#     u_prev = x[nx .+ (1:93)]
#     w = (u - u_prev) ./ h
#     J += 0.5 * 1.0 * dot(w[1:end-1], w[1:end-1])

# 	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal([1000.0 * ones(6); 100.0 * ones(nx-6)]) * (x[1:nx] - x_ref)
# 	J += 0.5 * transpose(u) * Diagonal([1.0 * ones(model.nu); zeros(nu - model.nu)]) * u
#     J += 1000.0 * u[end] # slack
#     J += 0.5 * transpose(u[model.nu + 8 .+ (1:40)]) * Diagonal(1.0e-3 * ones(40)) * u[model.nu + 8 .+ (1:40)]

# 	return J
# end

# function objT(x, u, w)
# 	J = 0.0
# 	J += 0.5 * transpose(x[1:nx] - xT) * Diagonal([1000.0 * ones(6); 1000.0 * ones(nx-6)]) * (x[1:nx] - xT)
#     J -= 100.0 * x[36 + model.nu + 5]
#     return J
# end

# c1 = DTO.Cost(obj1, nx, nu)
# ct = DTO.Cost(objt, nx + nθ , nu)
# cT = DTO.Cost(objT, nx + nθ , 0)
# obj = [c1, [ct for t = 2:T-1]..., cT];

# ## constraints
ql = q1
qu = q1

# initial condition
xl1 = [q1; q1]
xu1 = [q1; q1]
xlt = [-Inf * ones(nx); -Inf * ones(nθ)]
xut = [Inf * ones(nx); Inf * ones(nθ)]

# final condition
xlT = [-Inf * ones(nx); -Inf * ones(nθ)]
xuT = [Inf * ones(nx); Inf * ones(nθ)]

ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

bnd1 = DTO.Bound(nx, nu, state_lower=xl1, state_upper=xu1, action_lower=ul, action_upper=uu)
bndt = DTO.Bound(nx + nθ , nu, state_lower=xlt, state_upper=xut, action_lower=ul, action_upper=uu)
bndT = DTO.Bound(nx + nθ , 0, state_lower=xlT, state_upper=xuT)
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];

contact_constraints_equality(model, env, h, rand(nx), rand(nu), zeros(0))
contact_constraints_inequality_1(model, env, h, rand(nx), rand(nu), zeros(0))
contact_constraints_inequality_t(model, env, h, rand(nx + nθ ), rand(nu), zeros(0))
contact_constraints_inequality_T(model, env, h, rand(nx + nθ ), rand(nu), zeros(0))

function constraints_1(x, u, w)
    [
     # equality (32)
     contact_constraints_equality(model, env, h, x, u, w);
     # inequality (56)
     contact_constraints_inequality_1(model, env, h, x, u, w);
     x[6 .+ (1:12)] - q1[6 .+ (1:12)];
     x[18 + 6 .+ (1:12)] - q1[6 .+ (1:12)];
    ]
end

function constraints_t(x, u, w)
    [
     # equality (32)
     contact_constraints_equality(model, env, h, x, u, w);
     # inequality (64)
     contact_constraints_inequality_t(model, env, h, x, u, w);
     x[6 .+ (1:3)] - q1[6 .+ (1:3)];
     x[6 .+ (7:12)] - q1[6 .+ (7:12)];
     x[18 + 6 .+ (1:3)] - q1[6 .+ (1:3)];
     x[18 + 6 .+ (7:12)] - q1[6 .+ (7:12)];
    #  x[18 + 6 .+ (4:12)] - q1[6 .+ (4:12)];
    #  x[18 + 3] - q1[3];
    ]
end

function constraints_T(x, u, w)
    [
     # inequality (16)
     contact_constraints_inequality_T(model, env, h, x, u, w);
     x[6 .+ (1:3)] - qT[6 .+ (1:3)];
     x[6 .+ (7:12)] - qT[6 .+ (7:12)];

     x[18 .+ (1:3)] - qT[1:3];
     x[18 + 6 .+ (1:3)] - qT[6 .+ (1:3)];
     x[18 + 6 .+ (7:12)] - qT[6 .+ (7:12)];
     x[18 + 9 .+ 1] - qT[9 + 1];
     x[18 + 9 .+ 2] - qT[9 + 2];
     x[18 + 9 .+ 3] - qT[9 + 3];
     x[0 .+ (3:6)] - qT[3:6];
     x[18 .+ (3:6)] - qT[3:6];
    ]
end

con1 = DTO.Constraint(constraints_1, nx, nu, indices_inequality=collect(32 .+ (1:56)))
cont = DTO.Constraint(constraints_t, nx + nθ , nu, indices_inequality=collect(32 .+ (1:64)))
conT = DTO.Constraint(constraints_T, nx + nθ , nu, indices_inequality=collect(0 .+ (1:16)))
cons = [con1, [cont for t = 2:T-1]..., conT];

# ## problem
direct_solver = DTO.Solver(dyn, obj, cons, bnds,
    options=DTO.Options(
        tol=1.0e-3,
        constr_viol_tol=1.0e-3,
        max_iter=3000,
        ))

# ## initialize
x_interpolation = [x1, [[x1; zeros(nθ)] for t = 2:T]...]
u_guess = [1.0e-4 * rand(nu) for t = 1:T-1] # may need to run more than once to get good trajectory
DTO.initialize_states!(direct_solver, x_interpolation)
DTO.initialize_controls!(direct_solver, u_guess)

# ## solve
@time DTO.solve!(direct_solver)

# ## solution
x_sol, u_sol = DTO.get_trajectory(direct_solver)
@show x_sol[1]
@show x_sol[T]
maximum([u[nu] for u in u_sol[1:end-1]])

plot([x_sol[t][36 + model.nu + 5] for t = 2:T])
plot([x_sol[t][36 + model.nu + 6] for t = 2:T])

# plot([x[18 + 7] for x in x_sol], [x[18 + 8] for x in x_sol])
x_sol[end][36 + model.nu + 5]
x_sol[end][36 + model.nu + 6]
x_sol[end][36 + model.nu .+ (1:4)]

# ## visualize
vis = Visualizer()
render(vis)
# q_sol = state_to_configuration([x[1:nx] for x in x_sol])
visualize!(vis, model, [x_sol[1][1:18], [x[18 .+ (1:18)] for x in x_sol]...], Δt=h)

N_first = 10
N_last = 20
q_opt = [[x_sol[1][model.nq .+ (1:model.nq)] for t = 1:N_first]..., x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]..., [x_sol[end][model.nq .+ (1:model.nq)] for t = 1:N_last]...]
v_opt = [[(x_sol[1][model.nq .+ (1:model.nq)] - x_sol[1][0 .+ (1:model.nq)]) ./ h for t = 1:N_first]..., [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]..., [(x_sol[end][model.nq .+ (1:model.nq)] - x_sol[end][0 .+ (1:model.nq)]) ./ h for t = 1:N_last]...]
u_opt = [[u_sol[1][1:model.nu] for t = 1:N_first]..., [u[1:model.nu] for u in u_sol]..., [u_sol[end][1:model.nu] for t = 1:N_last]...]
γ_opt = [[u_sol[1][model.nu .+ (1:8)] for t = 1:N_first]..., [u[model.nu .+ (1:8)] for u in u_sol]..., [u_sol[end][model.nu .+ (1:8)] for t = 1:N_last]...]
b_opt = [[u_sol[1][model.nu + 8 .+ (1:32)] for t = 1:N_first]..., [u[model.nu + 8 .+ (1:32)] for u in u_sol]..., [u_sol[end][model.nu + 8 .+ (1:32)] for t = 1:N_last]...]
ψ_opt = [[u_sol[1][model.nu + 8 + 32 .+ (1:8)] for t = 1:N_first]..., [u[model.nu + 8 + 32 .+ (1:8)] for u in u_sol]..., [u_sol[end][model.nu + 8 + 32 .+ (1:8)] for t = 1:N_last]...]
η_opt = [[u_sol[1][model.nu + 8 + 32 + 8 .+ (1:32)] for t = 1:N_first]..., [u[model.nu + 8 + 32 + 8 .+ (1:32)] for u in u_sol]..., [u_sol[end][model.nu + 8 + 32 + 8 .+ (1:32)] for t = 1:N_last]...]

qm = q_opt
vm = v_opt
um = u_opt
γm = γ_opt
bm = b_opt
ψm = ψ_opt
ηm = η_opt

μm = model.μ_world
hm = h
timesteps = range(0.0, stop=(h * (length(qm) - 2)), length=(length(qm) - 2))
plot(hcat(qm...)', labels="")
plot(timesteps, hcat(um...)', labels="")
plot(timesteps, hcat(γm...)', labels="")
plot(timesteps, hcat(bm...)', labels="")
plot(timesteps, hcat(ψm...)', labels="")
plot(timesteps, hcat(ηm...)', labels="")

using JLD2

@save joinpath(@__DIR__, "wall_stand_FL_4.jld2") qm um γm bm ψm ηm μm hm
@load joinpath(@__DIR__, "wall_stand_FL_4.jld2") qm um γm bm ψm ηm μm hm



# q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
# v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
# u_opt = [u[1:model.nu] for u in u_sol]
# λ_opt = [u[model.nu .+ (1:5)] for u in u_sol]
# b_opt = [u[model.nu + 5 .+ (1:20)] for u in u_sol]


# q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
# v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
# u_opt = [u[1:model.nu] for u in u_sol]
# γ_opt = [u[model.nu .+ (1:5)] for u in u_sol]
# b_opt = [u[model.nu + 5 .+ (1:20)] for u in u_sol]
# ψ_opt = [u[model.nu + 5 + 20 .+ (1:5)] for u in u_sol]
# η_opt = [u[model.nu + 5 + 20 + 5 .+ (1:20)] for u in u_sol]

# qm = q_opt
# vm = v_opt
# um = u_opt
# γm = γ_opt
# bm = b_opt
# ψm = ψ_opt
# ηm = η_opt

# μm = model.μ_world
# hm = h
# timesteps = range(0.0, stop=(h * (length(qm) - 2)), length=(length(qm) - 2))
# plot(hcat(qm...)', labels="")
# plot(timesteps, hcat(um...)', labels="")
# plot(timesteps, hcat(γm...)', labels="")
# plot(timesteps, hcat(bm...)', labels="")
# plot(timesteps, hcat(ψm...)', labels="")
# plot(timesteps, hcat(ηm...)', labels="")

# using JLD2
# @save joinpath(@__DIR__, "wall_stand_FL.jld2") qm um γm bm ψm ηm μm hm
# @load joinpath(@__DIR__, "wall_stand_FL.jld2") qm um γm bm ψm ηm μm hm
