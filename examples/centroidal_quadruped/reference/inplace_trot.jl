# ## model 
include("ci_model.jl") 

# ## horizon 

freq = 100
h = 1.0 / freq
T = Int(floor(0.5 * 0.65 / h)) + 1
Tm = 17 # T / 2

# ## centroidal_quadruped 
# https://github.com/ChiyenLee/QuadrupedBalance.jl/blob/dojo/notebooks/Centroidal_model.jl.ipynb
model = RoboDojo.centroidal_quadruped
model.inertia_body = Matrix(Diagonal([0.05, 0.25, 0.3]))
model.mass_body = 12.54 # mass
model.mass_foot = 12.54 / 100

nx = 2 * model.nq
nc = 4 #model.nc
nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
nθ = 5 

RoboDojo.mass_matrix(model, ones(model.nq))
RoboDojo.dynamics_bias(model, ones(model.nq), ones(model.nq))
RoboDojo.contact_jacobian(model, ones(model.nq))[1:12, :]

# ## model
d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dyn1(model, RoboDojo.mass_matrix, RoboDojo.dynamics_bias, [h], y, x, u, w), nx + nθ + nx + model.nu, nx, nu)
dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(model, RoboDojo.mass_matrix, RoboDojo.dynamics_bias, [h], y, x, u, w), nx + nθ + nx + model.nu, nx + nθ + nx + model.nu, nu)

dyn = [d1, [dt for t = 2:T-1]...]

# ## initial conditions
mode = :right
foot_height = 0.05
q1 = nominal_configuration(model) 
qM = nominal_configuration(model)
if mode == :left
    # qM[6 + 3] += foot_height # front left
    # qM[15 + 3] += foot_height # back right
elseif mode == :right
    qM[9 + 3] += foot_height # front right
    qM[12 + 3] += foot_height # back left
end

qT = nominal_configuration(model)

RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, [qM], Δt=h);

q_ref = [linear_interpolation(q1, qM, Tm)..., linear_interpolation(qM, qT, Tm)...]
x_ref = [[q_ref[t]; q_ref[t+1]] for t = 1:T]
x1 = x_ref[1]
xM = x_ref[Tm]
xT = x_ref[T]

# RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, q_ref, Δt=h);

# ## objective
obj = DTO.Cost{Float64}[] 
for t = 1:T 
    if t == T 
        function objT(x, u, w)
            J = 0.0 
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-1 * dot(v, v)
            J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t]) 
            return J 
        end
        push!(obj, DTO.Cost(objT, nx + nθ + nx + model.nu, 0))
    elseif t == 1 
        function obj1(x, u, w)
            J = 0.0 
            v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
            J += 0.5 * 1.0e-3 * dot(v, v)
            J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t]) 
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
            u_previous = x[nx + 5 + nx .+ (1:model.nu)] 
            u_control = u[1:model.nu] 
            w = (u_control - u_previous) ./ h
            J += 0.5 * 1.0e-2 * dot(w, w)
            J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(1000.0 * ones(nx)) * (x[1:nx] - x_ref[t]) 
            J += 0.5 * transpose(u[1:model.nu]) * Diagonal(1.0e-3 * ones(model.nu)) * u[1:model.nu]
            J += 1000.0 * u[end] # slack
            return J 
        end
        push!(obj, DTO.Cost(objt, nx + nθ + nx + model.nu, nu))
    end 
end

# ## constraints
ql = q1
qu = q1

# initial condition
xl1 = [q1; q1] 
xu1 = [q1; q1]
xlt = [-Inf * ones(nx); -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)] 
xut = [Inf * ones(nx); Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

# final condition
xlT = [qT; qT; -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)] 
xuT = [qT; qT; Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
bndt = DTO.Bound(nx + nθ + nx + model.nu, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
bndT = DTO.Bound(nx + nθ + nx + model.nu, 0, xl=xlT, xu=xuT)
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];

cons = DTO.Constraint{Float64}[]
for t = 1:T     
    if t == 1
        function constraints_1(x, u, w) 
            [
            # equality (16)
            contact_constraints_equality(model, h, x, u, w); 
            # inequality (28)
            contact_constraints_inequality_1(model, h, x, u, w);

            # body/feet constraints
            # x[3] - x_ref[t][3]; # body height 
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            ]
        end
        push!(cons, DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(16 .+ (1:28))))
    elseif t == T 
        function constraints_T(x, u, w) 
            [
            # inequality (8)
            contact_constraints_inequality_T(model, h, x, u, w);

            # body/feet constraints
            # x[3] - x_ref[t][3]; # body height 
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            ]
        end
        push!(cons, DTO.Constraint(constraints_T, nx + nθ + nx + model.nu, nu, idx_ineq=collect(0 .+ (1:8))))
    else 
        function constraints_t(x, u, w) 
            [
            # equality (16)
            contact_constraints_equality(model, h, x, u, w); 
            # inequality (32)
            contact_constraints_inequality_t(model, h, x, u, w);

            # body/feet constraints
            # x[3] - x_ref[t][3]; # body height 
            # x[model.nq + 3] - x_ref[t][model.nq + 3]; # body height
            ]
        end
        push!(cons, DTO.Constraint(constraints_t, nx + nθ + nx + model.nu, nu, idx_ineq=collect(16 .+ (1:32))) )
    end
end

# ## problem 
tolerance = 1.0e-2
p = DTO.solver(dyn, obj, cons, bnds, 
    options=DTO.Options(
        tol=tolerance,
        constr_viol_tol=tolerance,
        ))

# ## initialize
x_interpolation = [x_ref[1], [[x_ref[t]; zeros(nθ); zeros(nx); zeros(model.nu)] for t = 2:T]...]
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
RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, x_sol, Δt=h);

q_opt = [x_sol[1][1:model.nq], [x[model.nq .+ (1:model.nq)] for x in x_sol]...]
v_opt = [(x[model.nq .+ (1:model.nq)] - x[0 .+ (1:model.nq)]) ./ h for x in x_sol]
u_opt = [u[1:model.nu] for u in u_sol] 
λ_opt = [u[model.nu .+ (1:4)] for u in u_sol]
b_opt = [u[model.nu + 4 .+ (1:16)] for u in u_sol]

# mirror 
q_mirror = Vector{Float64}[]
u_mirror = Vector{Float64}[] 
reflection = Diagonal([1.0; -1.0; 1.0]) 

for (t, q) in enumerate(q_opt)
    qfl = q[6 .+ (1:3)]  # front left
    qfr = q[9 .+ (1:3)]  # front right
    qbl = q[12 .+ (1:3)] # back left
    qbr = q[15 .+ (1:3)] # back right

    push!(q_mirror, [q[1:6]; reflection * qfr; reflection * qfl; reflection * qbr; reflection * qbl])

    if t < T 
        # set left using right 
        ufl = u_opt[t][0 .+ (1:3)]  # front left
        ufr = u_opt[t][3 .+ (1:3)] # front right
        ubl = u_opt[t][6 .+ (1:3)] # back left
        ubr = u_opt[t][9 .+ (1:3)] # back right
        if t < T 
            push!(u_mirror, [reflection * ufr; reflection * ufl; reflection * ubr; reflection * ubl])
        end
    end
end

using Plots
plot(hcat(q_opt..., q_mirror...)', xlabel="time", ylabel="configuration", label="") 
plot(hcat(u_opt..., u_mirror...)', xlabel="time", ylabel="control", label="") 

RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, [q_opt..., q_mirror...], Δt=h);
