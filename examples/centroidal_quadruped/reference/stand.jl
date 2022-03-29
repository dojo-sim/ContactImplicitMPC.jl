using RoboDojo 
using LinearAlgebra 
using DirectTrajectoryOptimization 
const DTO = DirectTrajectoryOptimization

function centroidal_quadruped_dyn(model, mass_matrix, dynamics_bias, h, y, x, u, w) 

    # dimensions
    nq = model.nq
    nu = model.nu 

    # configurations
    
    q1⁻ = x[1:nq] 
    q2⁻ = x[nq .+ (1:nq)]
    q2⁺ = y[1:nq]
    q3⁺ = y[nq .+ (1:nq)]

    # control 
    u_control = u[1:nu] 
    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:16)] 
    
    E = [1.0 0.0 -1.0 0.0; 
         0.0 1.0 0.0 -1.0] # friction mapping 
    J = RoboDojo.contact_jacobian(model, q2⁺)
    λ = transpose(J[1:12, :]) * [
                        [E * β[0  .+ (1:4)]; γ[1]];
                        [E * β[4  .+ (1:4)]; γ[2]];
                        [E * β[8  .+ (1:4)]; γ[3]];
                        [E * β[12 .+ (1:4)]; γ[4]]
                       ]
    [
     q2⁺ - q2⁻;
     RoboDojo.dynamics(model, q -> mass_matrix(model, q), (q, q̇) -> dynamics_bias(model, q, q̇), 
        h, q1⁻, q2⁺, u_control, zeros(model.nw), λ, q3⁺)
    ]
end

function centroidal_quadruped_dyn1(model, mass_matrix, dynamics_bias, h, y, x, u, w)
    nx = 2 * model.nq
    [
     centroidal_quadruped_dyn(model, mass_matrix, dynamics_bias, h, y, x, u, w);
     y[nx .+ (1:5)] - [u[model.nu .+ (1:4)]; u[end]];
    ]
end

function centroidal_quadruped_dynt(model, mass_matrix, dynamics_bias, h, y, x, u, w)
    nx = 2 * model.nq
    [
     centroidal_quadruped_dyn(model, mass_matrix, dynamics_bias, h, y, x, u, w);
     y[nx .+ (1:5)] - [u[model.nu .+ (1:4)]; u[end]];
    ]
end

function contact_constraints_inequality_1(model, h, x, u, w) 
    nq = model.nq
    nu = model.nu 
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    u_control = u[1:nu] 
    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:16)] 
    ψ = u[nu + 4 + 16 .+ (1:4)] 
    η = u[nu + 4 + 16 + 4 .+ (1:16)] 
    sα = u[nu + 4 + 16 + 4 + 16 .+ (1:1)]

    ϕ = RoboDojo.signed_distance(model, q3)[1:4]
  
    μ = model.friction_foot_world
    fc = μ .* γ[1:4] - [sum(β[0 .+ (1:4)]); sum(β[4 .+ (1:4)]); sum(β[8 .+ (1:4)]); sum(β[12 .+ (1:4)]);]

    [
     -ϕ; 
     -fc;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end

function contact_constraints_inequality_t(model, h, x, u, w) 
    nq = model.nq
    nu = model.nu 
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    u_control = u[1:nu] 
    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:16)] 
    ψ = u[nu + 4 + 16 .+ (1:4)] 
    η = u[nu + 4 + 16 + 4 .+ (1:16)] 
    sα = u[nu + 4 + 16 + 4 + 16 .+ (1:1)]

    ϕ = RoboDojo.signed_distance(model, q3)[1:4]
    γ⁻ = x[nx .+ (1:4)] 
    sα⁻ = x[nx + 4 .+ (1:1)]
    
    μ = model.friction_foot_world
    fc = μ .* γ[1:4] - [sum(β[0 .+ (1:4)]); sum(β[4 .+ (1:4)]); sum(β[8 .+ (1:4)]); sum(β[12 .+ (1:4)]);]

    [
     -ϕ; 
     -fc;
     γ⁻ .* ϕ .- sα⁻;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end


function contact_constraints_inequality_T(model, h, x, u, w) 
    nq = model.nq
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    ϕ = RoboDojo.signed_distance(model, q3)[1:4]
    γ⁻ = x[nx .+ (1:4)] 
    sα⁻ = x[nx + 4 .+ (1:1)]
   
    [
     -ϕ; 
     γ⁻ .* ϕ .- sα⁻;
    ]
end

function contact_constraints_equality(model, h, x, u, w) 
    nq = model.nq
    nu = model.nu 

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:16)] 
    ψ = u[nu + 4 + 16 .+ (1:4)] 
    η = u[nu + 4 + 16 + 4 .+ (1:16)] 
    sα = u[nu + 4 + 16 + 4 + 16 .+ (1:1)]
   
    E = [1.0 0.0 -1.0 0.0; 
         0.0 1.0 0.0 -1.0]
    v = (q3 - q2) ./ h[1]
    vT = vcat([E' * v[(i-1) * 3 .+ (1:2)] for i = 1:4]...)
    
    ψ_stack = vcat([ψi * ones(4) for ψi in ψ]...)
    
    [
     η - vT - ψ_stack;
    ]
end

# ## horizon 
T = 11 
h = 0.1

# ## centroidal_quadruped 
model = RoboDojo.centroidal_quadruped
nx = 2 * model.nq
nc = 4 #model.nc
nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
nθ = 5 

RoboDojo.mass_matrix(model, ones(model.nq))
RoboDojo.dynamics_bias(model, ones(model.nq), ones(model.nq))
RoboDojo.contact_jacobian(model, ones(model.nq))[1:12, :]

# ## model
d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dyn1(model, RoboDojo.mass_matrix, RoboDojo.dynamics_bias, [h], y, x, u, w), nx + nθ, nx, nu)
dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(model, RoboDojo.mass_matrix, RoboDojo.dynamics_bias, [h], y, x, u, w), nx + nθ, nx + nθ, nu)

dyn = [d1, [dt for t = 2:T-1]...]

# ## initial conditions
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
	J += 0.5 * transpose(x - x_ref) * Diagonal(ones(nx)) * (x - x_ref) 
	J += 0.5 * transpose(u) * Diagonal([ones(model.nu); zeros(nu - model.nu)]) * u
    J += 1000.0 * u[end] # slack
	return J
end

function objt(x, u, w)
	J = 0.0 
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal(ones(nx)) * (x[1:nx] - x_ref) 
	J += 0.5 * transpose(u) * Diagonal([ones(model.nu); zeros(nu - model.nu)]) * u
    J += 1000.0 * u[end] # slack
	return J
end

function objT(x, u, w)
	J = 0.0 
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal(ones(nx)) * (x[1:nx] - x_ref) 
    return J
end

c1 = DTO.Cost(obj1, nx, nu)
ct = DTO.Cost(objt, nx + nθ, nu)
cT = DTO.Cost(objT, nx + nθ, 0)
obj = [c1, [ct for t = 2:T-1]..., cT];

# ## constraints

ql = q1
qu = q1

# initial condition
xl1 = [q1; q1] 
xu1 = [q1; q1]
xlt = [-Inf * ones(nx); -Inf * ones(nθ)] 
xut = [Inf * ones(nx); Inf * ones(nθ)]

# final condition
xlT = [q1; q1; -Inf * ones(nθ)] 
xuT = [q1; q1; Inf * ones(nθ)]

ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
bndt = DTO.Bound(nx + nθ, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
bndT = DTO.Bound(nx + nθ, 0, xl=xlT, xu=xuT)
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];


function constraints_1(x, u, w) 
    [
     # equality (16)
     contact_constraints_equality(model, h, x, u, w); 
     # inequality (28)
     contact_constraints_inequality_1(model, h, x, u, w);
    ]
end

function constraints_t(x, u, w) 
    [
     # equality (16)
     contact_constraints_equality(model, h, x, u, w); 
     # inequality (32)
     contact_constraints_inequality_t(model, h, x, u, w);
    ]
end

function constraints_T(x, u, w) 
    [
     # inequality (8)
     contact_constraints_inequality_T(model, h, x, u, w);
    ]
end

con1 = DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(16 .+ (1:28))) 
cont = DTO.Constraint(constraints_t, nx + nθ, nu, idx_ineq=collect(16 .+ (1:32))) 
conT = DTO.Constraint(constraints_T, nx + nθ, nu, idx_ineq=collect(0 .+ (1:8))) 
cons = [con1, [cont for t = 2:T-1]..., conT];

# ## problem 
p = DTO.solver(dyn, obj, cons, bnds, 
    options=DTO.Options(
        tol=1.0e-2,
        constr_viol_tol=1.0e-2,
        ))

# ## initialize
x_interpolation = [x1, [[x1; zeros(nθ)] for t = 2:T]...]
u_guess = [1.0e-4 * rand(nu) for t = 1:T-1] # may need to run more than once to get good trajectory
DTO.initialize_states!(p, x_interpolation)
DTO.initialize_controls!(p, u_guess)

# ## solve
@time DTO.solve!(p)

# ## solution
x_sol, u_sol = DTO.get_trajectory(p)
@show x_sol[1]
@show x_sol[T]
sum([u[nu] for u in u_sol[1:end-1]])

# ## visualize 
vis = Visualizer() 
render(vis)
# q_sol = state_to_configuration([x[1:nx] for x in x_sol])
RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, x_sol, Δt=h);
