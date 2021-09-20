ENV["GUROBI_HOME"] = "/home/simon/software/gurobi912/linux64/"
ENV["GRB_LICENSE_FILE"] = "/home/simon/software/gurobi912/gurobi.lic"

using Gurobi
using Convex
using Random

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 100)
set_optimizer_attribute(model, "Presolve", 0)



using Convex
using LinearAlgebra
using SparseArrays
using Test
using Plots

##
# Data taken from http://people.sc.fsu.edu/~jburkardt/datasets/knapsack_01/knapsack_01.html
w = [23; 31; 29; 44; 53; 38; 63; 85; 89; 82]
C = 165
p =  [92; 57; 49; 68; 60; 43; 67; 84; 87; 72];
n = length(w)

# u = Variable(2)
x = Variable(n, :Bin)
problem = maximize(dot(p, x), dot(w, x) <= C)
solve!(problem, Gurobi.Optimizer)
evaluate(x)


##
aux(str) = joinpath(@__DIR__, "aux_files", str) # path to auxiliary files
include(aux("antidiag.jl"))

n = 8
x = Variable((n, n), :Bin)
# At most one queen on any anti-diagonal
constr = Constraint[sum(antidiag(x, k)) <= 1 for k = -n+2:n-2]
# At most one queen on any diagonal
constr += Constraint[sum(diag(x, k)) <= 1 for k = -n+2:n-2]
# Exactly one queen per row and one queen per column
constr += Constraint[sum(x, dims=1) == 1, sum(x, dims=2) == 1]
p = satisfy(constr)
solve!(p, Gurobi.Optimizer)


for k = -n+2:n-2
	@test evaluate(sum(antidiag(x, k))) <= 1
	@test evaluate(sum(diag(x, k))) <= 1
end
@test all(evaluate(sum(x, dims=1)) .≈ 1)
@test all(evaluate(sum(x, dims=2)) .≈ 1)









ENV["GUROBI_HOME"] = "/home/simon/software/gurobi912/linux64/"
ENV["GRB_LICENSE_FILE"] = "/home/simon/software/gurobi912/gurobi.lic"

using Gurobi
using JuMP
using LinearAlgebra
using Random
using Plots

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 100)
set_optimizer_attribute(model, "Presolve", 0)


##
# MILP
Random.seed!(314)

# number of clients
m = 12
# number of facility locations
n = 5

# Clients' locations
Xc = rand(m)
Yc = rand(m)

# Facilities' potential locations
Xf = rand(n)
Yf = rand(n)

# Fixed costs
f = ones(n);

# Distance
c = zeros(m, n)
for i in 1:m
    for j in 1:n
        c[i, j] = LinearAlgebra.norm([Xc[i] - Xf[j], Yc[i] - Yf[j]], 2)
    end
end

Plots.scatter(
    Xc,
    Yc,
    label = "Clients",
    markershape = :circle,
    markercolor = :blue,
)
Plots.scatter!(
    Xf,
    Yf,
    label = "Facility",
    markershape = :square,
    markercolor = :white,
    markersize = 6,
    markerstrokecolor = :red,
    markerstrokewidth = 2,
)


ufl = Model(Gurobi.Optimizer)


@variable(ufl, y[1:n], Bin);
@variable(ufl, v[1:n]);
@variable(ufl, x[1:m, 1:n], Bin);

@constraint(ufl, client_service[i in 1:m], sum(x[i, j] for j in 1:n) == 1);
@constraint(ufl, open_facility[i in 1:m, j in 1:n], x[i, j] <= y[j]);
@objective(ufl, Min, sum(v.^2) + f'y + sum(c .* x));

optimize!(ufl)
println("Optimal value: ", objective_value(ufl))

c .* x

x_ = value.(x) .> 1 - 1e-5
y_ = value.(y) .> 1 - 1e-5

p = Plots.scatter(
    Xc,
    Yc,
    markershape = :circle,
    markercolor = :blue,
    label = nothing,
)

mc = [(y_[j] ? :red : :white) for j in 1:n]
Plots.scatter!(
    Xf,
    Yf,
    markershape = :square,
    markercolor = mc,
    markersize = 6,
    markerstrokecolor = :red,
    markerstrokewidth = 2,
    label = nothing,
)


for i in 1:m
    for j in 1:n
        if x_[i, j] == 1
            Plots.plot!(
                [Xc[i], Xf[j]],
                [Yc[i], Yf[j]],
                color = :black,
                label = nothing,
            )
        end
    end
end

p




#
# ##
# #MIQP
#
# T = 1
# N = 3
# n = 4
# m = 3
# A = [Matrix(i * Diagonal(ones(n))) for i = 1:N]
# B = [Matrix(i * Diagonal(ones(m))) for i = 1:N]
#
# P = Matrix(Diagonal(ones(n)))
# Q = Matrix(Diagonal(ones(n)))
# R = Matrix(Diagonal(ones(m)))
# x0 = ones(n)
#
# b = Variable(N, :Bin)
# x = Variable(T+1)
# u = Variable(T)
#
# obj = sum(b)
# con = sum(b) >= 1
# prob = minimize(obj, con)
# solve!(prob, GLPK.Optimizer)





##
# MIQP
Random.seed!(314)
T = 1
N = 3
n = 4
m = 3
A = [Matrix(i * Diagonal(ones(n))) for i = 1:N]
B = [rand(n, m) for i = 1:N]

P = Matrix(Diagonal(ones(n)))
Q = Matrix(Diagonal(ones(n)))
R = Matrix(Diagonal(ones(m)))
x0 = 10*ones(n)

ufl = Model(Gurobi.Optimizer)

@variable(ufl, b[1:N], Bin);
@variable(ufl, x[1:T+1, 1:n]);
@variable(ufl, u[1:m]);

@constraint(ufl, initial_state, x[1,:] .== x0);
@constraint(ufl, dynamics[t in 1:T], x[t+1,:] .== A[1] * x[t,:] + B[1] * u);
@objective(ufl, Min, sum(x.^2) + sum(u.^2));

optimize!(ufl)
println("Optimal value: ", objective_value(ufl))

x_ = value.(x)
u_ = value.(u)


obj = sum(b)
con = sum(b) >= 1
prob = minimize(obj, con)
solve!(prob, GLPK.Optimizer)



# number of clients
m = 12
# number of facility locations
n = 5

# Clients' locations
Xc = rand(m)
Yc = rand(m)

# Facilities' potential locations
Xf = rand(n)
Yf = rand(n)

# Fixed costs
f = ones(n);

# Distance
c = zeros(m, n)
for i in 1:m
    for j in 1:n
        c[i, j] = LinearAlgebra.norm([Xc[i] - Xf[j], Yc[i] - Yf[j]], 2)
    end
end

Plots.scatter(
    Xc,
    Yc,
    label = "Clients",
    markershape = :circle,
    markercolor = :blue,
)
Plots.scatter!(
    Xf,
    Yf,
    label = "Facility",
    markershape = :square,
    markercolor = :white,
    markersize = 6,
    markerstrokecolor = :red,
    markerstrokewidth = 2,
)


ufl = Model(Gurobi.Optimizer)


@variable(ufl, y[1:n], Bin);
@variable(ufl, v[1:n]);
@variable(ufl, x[1:m, 1:n], Bin);

@constraint(ufl, client_service[i in 1:m], sum(x[i, j] for j in 1:n) == 1);
@constraint(ufl, open_facility[i in 1:m, j in 1:n], x[i, j] <= y[j]);
@objective(ufl, Min, sum(v.^2) + f'y + sum(c .* x));

optimize!(ufl)
println("Optimal value: ", objective_value(ufl))

c .* x

x_ = value.(x) .> 1 - 1e-5
y_ = value.(y) .> 1 - 1e-5

p = Plots.scatter(
    Xc,
    Yc,
    markershape = :circle,
    markercolor = :blue,
    label = nothing,
)

mc = [(y_[j] ? :red : :white) for j in 1:n]
Plots.scatter!(
    Xf,
    Yf,
    markershape = :square,
    markercolor = mc,
    markersize = 6,
    markerstrokecolor = :red,
    markerstrokewidth = 2,
    label = nothing,
)


for i in 1:m
    for j in 1:n
        if x_[i, j] == 1
            Plots.plot!(
                [Xc[i], Xf[j]],
                [Yc[i], Yf[j]],
                color = :black,
                label = nothing,
            )
        end
    end
end

p



##

using Gurobi
using JuMP
using LinearAlgebra
using MeshCat
using Plots
using Test

# dimensions
T = 2 # horizon (T+1 states, T controls)
n = 2 # state dim
m = 1 # control dim
nd = 2 # number of doamins

# dynamics parameters
mass = 1.0
l = 1.0
g = 10.0
@warn "changed"
k = 100.0
d = 0.1
t_s = 0.01

include("structures.jl")
include("visuals.jl")

# dynamics model
A1 = t_s * I + [0   1;
	  g/l 0]
B1 = reshape([0, 1 / (mass*l^2)], (n,m))
c1 = [0, 0]

A2 = t_s * I + [0         1;
	  g/l - k/mass 0]
B2 = reshape([0, 1 / (mass*l^2)], (n,m))
c2 = [0, k * d / (mass * l)]

A = [A1, A2, A1, A2]
B = [B1, B2, B1, B2]
c = [c1, c2, c1, c2]

# domains
x_min = [-d/l*2, -1.5]
x_max = [ d/l*2,  1.5]

x_min_1 = [-d/l*2, -1.5]
x_max_1 = [ d/l,    1.5]

x_min_2 = [ d/l,   -1.5]
x_max_2 = [ d/l*2,  1.5]

u_min = [-4.0]
u_max = [ 4.0]

C0 = BoxDomain12(x_min, x_max, u_min, u_max)
C1 = BoxDomain12(x_min_1, x_max_1, u_min, u_max)
C2 = BoxDomain12(x_min_2, x_max_2, u_min, u_max)
C3 = BoxDomain12([-1e5, -1e5], [d/l,  1e5], [-1e5], [1e5])
C4 = BoxDomain12([ d/l, -1e5], [1e5,  1e5], [-1e5], [1e5])
CS = [C1.S, C2.S]
CR = [C1.R, C2.R]
CT = [C1.T, C2.T]
C = [C1, C2, C3, C4]

# Extreme values
β = 1e4
Mstar = [β * ones(2n+2m) for i=1:nd]
MU = β * ones(n) # upper bound
ML = - β * ones(n) # lower bound

# Initial state
x0 = [0.02, 0.0]

# Cost function
Q = 0.0
R = 1.0

# MIQP problem
model = Model(Gurobi.Optimizer)

# variables
@variable(model, x[1:T+1, 1:n])
@variable(model, u[1:T, 1:m])
@variable(model, δ[1:nd, 1:T], Bin)
@variable(model, z[1:nd, 1:T, 1:n])

# constraints
# intitial state
@constraint(model, initial_state, x[1,:] .== x0);
# 16.a
@constraint(model, c16a[i in 1:nd, t in 1:T], CS[i] * x[t,:] + CR[i] * u[t,:] - CT[i] .<= Mstar[i] * (1 - δ[i,t]));
# 16.b
@constraint(model, c16b[t in 1:T], sum(δ[:, t]) == 1)
# 18
@constraint(model, c18[t in 1:T], x[t+1] .== sum([z[i,t,:] for i = 1:nd]))
# 22.a
@constraint(model, c22a[i in 1:nd, t in 1:T], z[i,t,:] .<= MU * δ[i,t])
# 22.b
@constraint(model, c22b[i in 1:nd, t in 1:T], z[i,t,:] .>= ML * δ[i,t])
# 22.c
@constraint(model, c22c[i in 1:nd, t in 1:T], z[i,t,:] .<= A[i] * x[t,:] + B[i] * u[t,:] + c[i] - ML * (1 - δ[i,t]))
# 22.d
@constraint(model, c22d[i in 1:nd, t in 1:T], z[i,t,:] .>= A[i] * x[t,:] + B[i] * u[t,:] + c[i] - MU * (1 - δ[i,t]))

# objective
# @objective(model, Min, Q * sum(x.^2) + R * sum(u.^2))
@objective(model, Min, R * sum(u.^2))

# solve
optimize!(model)
value.(x)
value.(u)
value.(δ)
value.(z)


println("Optimal value: ", objective_value(model))





H = 100
x0 = [0.0, 0.0]
u = [1e-6 * ones(m) for t = 1:H]
dt = 0.01

p = OpenLoopPolicy16(u, dt = dt)
sim = Simulator12(p, dt = dt)

simulate!(sim, x0, H)
DTx = range(0.0, step = dt, length = size(sim.x)[1])
DTu = range(0.0, step = dt, length = size(sim.u)[1])
plot(DTx, hcat(sim.x...)', label = "state")
plot(DTu, hcat(sim.u...)', label = "control")


# vis = Visualizer()
# open(vis)
build_robot!(vis, nothing)
set_robot!(vis, nothing, [0.1, 0.0])

function kinematics(x::AbstractVector{T}) where {T}
	θ = x[1] + π/2
	θd = x[2]
	p = l * [cos(θ), sin(θ)]
	return p
end

visualize_robot!(vis, pushbot, sim.x)
