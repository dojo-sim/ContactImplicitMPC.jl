using RoboDojo 
using LinearAlgebra 
using DirectTrajectoryOptimization 
const DTO = DirectTrajectoryOptimization

function hopper_dyn(mass_matrix, dynamics_bias, h, y, x, u, w) 
    model = RoboDojo.hopper

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
    β = u[nu + 4 .+ (1:4)] 
    
    E = [1.0 -1.0] # friction mapping 
    J = RoboDojo.contact_jacobian(model, q2⁺)
    λ = transpose(J) * [[E * β[1:2]; γ[1]];
                        [E * β[3:4]; γ[2]];
                         γ[3:4]]
    λ[3] += (model.body_radius * E * β[1:2])[1] # friction on body creates a moment

    [
     q2⁺ - q2⁻;
     RoboDojo.dynamics(model, mass_matrix, dynamics_bias, 
        h, q1⁻, q2⁺, u_control, zeros(model.nw), λ, q3⁺)
    ]
end

function hopper_dyn1(mass_matrix, dynamics_bias, h, y, x, u, w)
    nx = 8 
    [
     hopper_dyn(mass_matrix, dynamics_bias, h, y, x, u, w);
     y[nx .+ (1:5)] - [u[2 .+ (1:4)]; u[end]];
     y[nx + 5 .+ (1:nx)] - x
    ]
end

function hopper_dynt(mass_matrix, dynamics_bias, h, y, x, u, w)
    nx = 8
    [
     hopper_dyn(mass_matrix, dynamics_bias, h, y, x, u, w);
     y[nx .+ (1:5)] - [u[2 .+ (1:4)]; u[end]];
     y[nx + 5 .+ (1:nx)] - x[nx + 5 .+ (1:nx)]
    ]
end

function contact_constraints_inequality_1(h, x, u, w) 
    model = RoboDojo.hopper

    nq = model.nq
    nu = model.nu 
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    u_control = u[1:nu] 
    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:4)] 
    ψ = u[nu + 4 + 4 .+ (1:2)] 
    η = u[nu + 4 + 4 + 2 .+ (1:4)] 
    sα = u[nu + 4 + 4 + 2 + 4 .+ (1:1)]

    ϕ = RoboDojo.signed_distance(model, q3) 
  
    μ = [model.friction_body_world; model.friction_foot_world]
    fc = μ .* γ[1:2] - [sum(β[1:2]); sum(β[3:4])]

    [
     -ϕ; 
     -fc;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end

function contact_constraints_inequality_t(h, x, u, w) 
    model = RoboDojo.hopper

    nq = model.nq
    nu = model.nu 
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:4)] 
    ψ = u[nu + 4 + 4 .+ (1:2)] 
    η = u[nu + 4 + 4 + 2 .+ (1:4)] 
    sα = u[nu + 4 + 4 + 2 + 4 .+ (1:1)]

    ϕ = RoboDojo.signed_distance(model, q3) 
    γ⁻ = x[nx .+ (1:4)] 
    sα⁻ = x[nx + 4 .+ (1:1)]
    
    μ = [model.friction_body_world; model.friction_foot_world]
    fc = μ .* γ[1:2] - [sum(β[1:2]); sum(β[3:4])]

    [
     -ϕ; 
     -fc;
     γ⁻ .* ϕ .- sα⁻;
     β .* η .- sα;
     ψ .* fc  .- sα;
    ]
end

function contact_constraints_inequality_T(h, x, u, w) 
    model = RoboDojo.hopper

    nq = model.nq
    nx = 2nq

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    ϕ = RoboDojo.signed_distance(model, q3) 
    γ⁻ = x[nx .+ (1:4)] 
    sα⁻ = x[nx + 4 .+ (1:1)]
   
    [
     -ϕ; 
     γ⁻ .* ϕ .- sα⁻;
    ]
end

function contact_constraints_equality(h, x, u, w) 
    model = RoboDojo.hopper

    nq = model.nq
    nu = model.nu 

    q2 = x[1:nq] 
    q3 = x[nq .+ (1:nq)] 

    γ = u[nu .+ (1:4)] 
    β = u[nu + 4 .+ (1:4)] 
    ψ = u[nu + 4 + 4 .+ (1:2)] 
    η = u[nu + 4 + 4 + 2 .+ (1:4)] 
   
    v = (q3 - q2) ./ h[1]
    vT_body = v[1] + model.body_radius * v[3]
    vT_foot = (RoboDojo.kinematics_foot_jacobian(model, q3) * v)[1]
    vT = [vT_body; -vT_body; vT_foot; -vT_foot]
    
    ψ_stack = [ψ[1] * ones(2); ψ[2] * ones(2)]
    
    [
     η - vT - ψ_stack;
    ]
end

# ## horizon 
T = 21 
h = 0.05

# ## hopper 
nx = 2 * RoboDojo.hopper.nq
nu = RoboDojo.hopper.nu + 4 + 4 + 2 + 4 + 1

# ## model
mass_matrix, dynamics_bias = RoboDojo.codegen_dynamics(RoboDojo.hopper)
d1 = DTO.Dynamics((y, x, u, w) -> hopper_dyn1(mass_matrix, dynamics_bias, [h], y, x, u, w), 2 * nx + 5, nx, nu)
dt = DTO.Dynamics((y, x, u, w) -> hopper_dynt(mass_matrix, dynamics_bias, [h], y, x, u, w), 2 * nx + 5, 2 * nx + 5, nu)

dyn = [d1, [dt for t = 2:T-1]...];

# ## initial conditions
q1 = [0.0; 0.5 + RoboDojo.hopper.foot_radius; 0.0; 0.5]
qM = [0.5; 0.5 + RoboDojo.hopper.foot_radius; 0.0; 0.5]
qT = [1.0; 0.5 + RoboDojo.hopper.foot_radius; 0.0; 0.5]
q_ref = [0.5; 0.75 + RoboDojo.hopper.foot_radius; 0.0; 0.25]

x1 = [q1; q1]
xM = [qM; qM]
xT = [qT; qT]
x_ref = [q_ref; q_ref]

# ## gate 
GAIT = 1 
GAIT = 2 
GAIT = 3

if GAIT == 1 
	r_cost = 1.0e-1 
	q_cost = 1.0e-1
elseif GAIT == 2 
	r_cost = 1.0
	q_cost = 1.0
elseif GAIT == 3 
	r_cost = 1.0e-3
	q_cost = 1.0e-1
end

# ## objective
function obj1(x, u, w)
	J = 0.0 
	J += 0.5 * transpose(x - x_ref) * Diagonal([1.0; 10.0; 1.0; 10.0; 1.0; 10.0; 1.0; 10.0]) * (x - x_ref) 
	J += 0.5 * transpose(u) * Diagonal([r_cost * ones(RoboDojo.hopper.nu); zeros(nu - RoboDojo.hopper.nu)]) * u
    J += 1000.0 * u[nu]
	return J
end

function objt(x, u, w)
	J = 0.0 
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal(q_cost * [1.0; 10.0; 1.0; 10.0; 1.0; 10.0; 1.0; 10.0]) * (x[1:nx] - x_ref)
	J += 0.5 * transpose(u) * Diagonal([r_cost * ones(RoboDojo.hopper.nu); zeros(nu - RoboDojo.hopper.nu)]) * u
    J += 1000.0 * u[nu]
	return J
end

function objT(x, u, w)
	J = 0.0 
	J += 0.5 * transpose(x[1:nx] - x_ref) * Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]) * (x[1:nx] - x_ref)
    return J
end

c1 = DTO.Cost(obj1, nx, nu)
ct = DTO.Cost(objt, 2 * nx + 5, nu)
cT = DTO.Cost(objT, 2 * nx + 5, 0)
obj = [c1, [ct for t = 2:T-1]..., cT];

# ## constraints

ql = [-Inf; 0; -Inf; 0.0]
qu = [Inf; Inf; Inf; 1.0]
xl1 = [q1; ql] 
xu1 = [q1; qu]
xlt = [ql; ql; -Inf * ones(5); -Inf * ones(nx)] 
xut = [qu; qu; Inf * ones(5); Inf * ones(nx)]
ul = [-10.0; -10.0; zeros(nu - 2)]
uu = [10.0; 10.0; Inf * ones(nu - 2)]

bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
bndt = DTO.Bound(2 * nx + 5, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
bndT = DTO.Bound(2 * nx + 5, 0, xl=xlt, xu=xut)
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];

function constraints_1(x, u, w) 
    [
     # equality (8)
   	 RoboDojo.kinematics_foot(RoboDojo.hopper, x[1:RoboDojo.hopper.nq]) - RoboDojo.kinematics_foot(RoboDojo.hopper, x1[1:RoboDojo.hopper.nq]);
	 RoboDojo.kinematics_foot(RoboDojo.hopper, x[RoboDojo.hopper.nq .+ (1:RoboDojo.hopper.nq)]) - RoboDojo.kinematics_foot(RoboDojo.hopper, x1[RoboDojo.hopper.nq .+ (1:RoboDojo.hopper.nq)]);
     contact_constraints_equality(h, x, u, w); 
     # inequality (12)
     contact_constraints_inequality_1(h, x, u, w);
    ]
end

function constraints_t(x, u, w) 
    [
     # equality (4)
     contact_constraints_equality(h, x, u, w); 
     # inequality (16)
     contact_constraints_inequality_t(h, x, u, w);
    ]
end

function constraints_T(x, u, w) 
    x_travel = 0.5
    θ = x[nx + 5 .+ (1:nx)]
    [
     # equality (6)
	 x[1:RoboDojo.hopper.nq][collect([2, 3, 4])] - θ[1:RoboDojo.hopper.nq][collect([2, 3, 4])];
	 x[RoboDojo.hopper.nq .+ (1:RoboDojo.hopper.nq)][collect([2, 3, 4])] - θ[RoboDojo.hopper.nq .+ (1:RoboDojo.hopper.nq)][collect([2, 3, 4])];
     # inequality (10)
     x_travel - (x[1] - θ[1])
	 x_travel - (x[RoboDojo.hopper.nq + 1] - θ[RoboDojo.hopper.nq + 1])
     contact_constraints_inequality_T(h, x, u, w);
    ]
end

con1 = DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(8 .+ (1:12))) 
cont = DTO.Constraint(constraints_t, 2nx + 5, nu, idx_ineq=collect(4 .+ (1:16))) 
conT = DTO.Constraint(constraints_T, 2nx + 5, nu, idx_ineq=collect(6 .+ (1:10))) 
cons = [con1, [cont for t = 2:T-1]..., conT];

# ## problem 
p = DTO.solver(dyn, obj, cons, bnds, 
    options=DTO.Options(
        tol=1.0e-2,
        constr_viol_tol=1.0e-2,
        print_level=0))

# ## initialize
x_interpolation = [x1, [[x1; zeros(5); x1] for t = 2:T]...]
u_guess = [[0.0; RoboDojo.hopper.gravity * RoboDojo.hopper.mass_body * 0.5 * h[1]; 1.0e-1 * ones(nu - 2)] for t = 1:T-1] # may need to run more than once to get good trajectory
DTO.initialize_states!(p, x_interpolation)
DTO.initialize_controls!(p, u_guess);

# ## solve
@time DTO.solve!(p)

# ## solution
x_sol, u_sol = DTO.get_trajectory(p)
@show x_sol[1]
@show x_sol[T]
sum([u[nu] for u in u_sol[1:end-1]])
x_sol[1] - x_sol[T][1:nx]

# ## visualize 
vis = Visualizer() 
render(vis)
# q_sol = state_to_configuration([x[1:nx] for x in x_sol])
RoboDojo.visualize!(vis, RoboDojo.hopper, x_sol, Δt=h);
