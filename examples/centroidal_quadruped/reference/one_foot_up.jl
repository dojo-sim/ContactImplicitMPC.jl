# ## model
include("trajopt_model.jl")

vis = Visualizer()
open(vis)

# ## horizon
T = 101
Tm = 51
h = 0.01

# ## centroidal_quadruped
s = get_simulation("centroidal_quadruped", "flat_3D_lc", "flat")
model = s.model
env = s.env
nx = 2 * model.nq
nc = 4 # model.nc
nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
nθ = 5

# ## model
d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dyn1(
	model, env, [h], y, x, u, w), nx + nθ + nx + model.nu, nx, nu)
dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(
	model, env, [h], y, x, u, w), nx + nθ + nx + model.nu, nx + nθ + nx + model.nu, nu)

dyn = [d1, [dt for t = 2:T-1]...]

# ## initial conditions
body_height = 0.3
foot_x = 0.17
foot_y = 0.15
function nominal_configuration(model::CentroidalQuadruped)
    [
        0.0; 0.0; body_height;
        0.0; 0.0; 0.0;
        foot_x; foot_y; 0.0;
        foot_x;-foot_y; 0.0;
       -foot_x; foot_y; 0.0;
       -foot_x;-foot_y; 0.0;
    ]
end
q1 = nominal_configuration(model)
qM = nominal_configuration(model)
qT = nominal_configuration(model)
q_ref = nominal_configuration(model)

x1 = [q1; q1]
xM = [qM; qM]
xT = [qT; qT]
x_ref = [q_ref; q_ref]

# ## objective
function obj1(x, u, w)
	J = 0.0
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal(ones(nx)) * (x[1:nx] - x_ref)
	J += 0.5 * transpose(u) * Diagonal([1.0e-2 * ones(model.nu); zeros(nu - model.nu)]) * u
    vf1 = (x[18 .+ (7:9)] - x[7:9]) ./ h 
    J += 1.0 * dot(vf1, vf1)
    J += 1000.0 * u[end] # slack
	return J
end

function objt(x, u, w)
	J = 0.0
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal(ones(nx)) * (x[1:nx] - x_ref)
	J += 0.5 * transpose(u) * Diagonal([1.0e-2 * ones(model.nu); zeros(nu - model.nu)]) * u
    vf1 = (x[18 .+ (7:9)] - x[7:9]) ./ h 
    J += 1.0 * dot(vf1, vf1)
    J += 1000.0 * u[end] # slack
	return J
end

function objT(x, u, w)
	J = 0.0
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal(ones(nx)) * (x[1:nx] - x_ref)
    vf1 = (x[18 .+ (7:9)] - x[7:9]) ./ h 
    J += 1.0 * dot(vf1, vf1)
    return J
end

c1 = DTO.Cost(obj1, nx, nu)
ct = DTO.Cost(objt, nx + nθ + nx + model.nu, nu)
cT = DTO.Cost(objT, nx + nθ + nx + model.nu, 0)
obj = [c1, [ct for t = 2:T-1]..., cT];

# ## constraints
ql = q1
qu = q1

# initial condition
xl1 = [q1; q1]
xu1 = [q1; q1]
xlt = [-Inf * ones(nx); -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)]
xut = [Inf * ones(nx); Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

# final condition
xlT = [q1; q1; -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)]
xuT = [q1; q1; Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
bndt = DTO.Bound(nx + nθ + nx + model.nu, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
bndT = DTO.Bound(nx + nθ + nx + model.nu, 0, xl=xlT, xu=xuT)
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT]

function constraints_1(x, u, w)
    [
     # equality (16)
     contact_constraints_equality(model, env, h, x, u, w);
     # inequality (28)
     contact_constraints_inequality_1(model, env, h, x, u, w);
     (x[18 .+ (10:18)] - q1[10:18]);
    ]
end

function constraints_t(x, u, w)
    [
     # equality (16)
     contact_constraints_equality(model, env, h, x, u, w);
     # inequality (32)
     contact_constraints_inequality_t(model, env, h, x, u, w);
     (x[18 .+ (10:18)] - q1[10:18]);
    ]
end

function constraints_M(x, u, w)
    [
     # equality (16)
     contact_constraints_equality(model, env, h, x, u, w);
     # inequality (32)
     contact_constraints_inequality_t(model, env, h, x, u, w);
     x[9] - 0.1;
     (x[18 .+ (10:18)] - q1[10:18]);
    ]
end

function constraints_T(x, u, w)
    [
     # inequality (8)
     contact_constraints_inequality_T(model, env, h, x, u, w);
     (x[18 .+ (10:18)] - q1[10:18]);
    ]
end

function constraints_TT(x, u, w)
    [
     # equality (16)
     contact_constraints_equality(model, env, h, x, u, w);
     # inequality (32)
     contact_constraints_inequality_t(model, env, h, x, u, w);
     x[1:36] - xT;
    ]
end

con1 = DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(16 .+ (1:28)))
cont = DTO.Constraint(constraints_t, nx + nθ + nx + model.nu, nu, idx_ineq=collect(16 .+ (1:32)))
conM = DTO.Constraint(constraints_M, nx + nθ + nx + model.nu, nu, idx_ineq=collect(16 .+ (1:32)))
conTT = DTO.Constraint(constraints_TT, nx + nθ + nx + model.nu, nu, idx_ineq=collect(16 .+ (1:32)))
conT = DTO.Constraint(constraints_T, nx + nθ + nx + model.nu, nu, idx_ineq=collect(0 .+ (1:8)))
cons = [con1, [t == Tm ? conM : (t > T - 10 ? conTT : cont) for t = 2:T-1]..., conT];

# ## problem
p = DTO.solver(dyn, obj, cons, bnds,
    options=DTO.Options(
        tol=1.0e-2,
        constr_viol_tol=1.0e-2,
        ))

# ## initialize
x_interpolation = [x1, [[x1; zeros(nθ); zeros(nx); zeros(model.nu)] for t = 2:T]...]
u_guess = [1.0e-4 * rand(nu) for t = 1:T-1] # may need to run more than once to get good trajectory
DTO.initialize_states!(p, x_interpolation)
DTO.initialize_controls!(p, u_guess)

# ## solve
@time DTO.solve!(p)

# ## solution
x_sol, u_sol = DTO.get_trajectory(p)
@show x_sol[1]
@show x_sol[T]
maximum([u[nu] for u in u_sol[1:end-1]])

# ## visualize
# q_sol = state_to_configuration([x[1:nx] for x in x_sol])
visualize!(vis, model, x_sol, Δt=h);


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
plot(timesteps, hcat(qm[2:end-1]...)', labels="")
plot(timesteps, hcat(um...)', labels="")
plot(timesteps, hcat(γm...)', labels="")
plot(timesteps, hcat(bm...)', labels="")
plot(timesteps, hcat(ψm...)', labels="")
plot(timesteps, hcat(ηm...)', labels="")

visualize!(vis, model, qm, Δt=h);

using JLD2
@save joinpath(@__DIR__, "one_foot_up.jld2") qm um γm bm ψm ηm μm hm
@load joinpath(@__DIR__, "one_foot_up.jld2") qm um γm bm ψm ηm μm hm
