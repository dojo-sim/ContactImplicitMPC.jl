include(joinpath(@__DIR__, "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
# render(vis)

# get hopper model
# model_sim = get_model("quadruped", surf="sinusoidal")
model = get_model("quadruped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("quadruped", "gait1", load_type=:split_traj_alt, model=model)


plot_surface!(vis, model.env, ylims=[-0.5, 0.3])
anim = visualize_meshrobot!(vis, model, ref_traj.q)
visualize_robot!(vis, model, ref_traj.q, anim=anim)

# Test visualizer
plot_surface!(vis, model.env, ylims=[-0.5, 0.3])
build_robot!(vis, model, name=:quad0f)
build_robot!(vis, model, name=:quad0b)
mvis = build_meshrobot!(vis, model, name=:shadow_1, α=0.5)

t = 20
q = ref_traj.q[t] + [0,0,pi/1,0,0,0,0,0,0,0,0]
set_robot!(vis, model, q, name=:quad0f, offset=0.00)
set_robot!(vis, model, q, name=:quad0b, offset=0.264)
set_meshrobot!(Fvis, mvis, model, q, name=:shadow_1)



vis2 = Visualizer()
open(vis2)

matx = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0))
maty = MeshPhongMaterial(color = RGBA(0.5, 0.5, 0.5, 1.0))
matz = MeshPhongMaterial(color = RGBA(0.8, 0.8, 0.8, 1.0))
objx = Rect(Vec(0, 0, 0),Vec(0.50, 0.05, 0.05))
objy = Rect(Vec(0, 0, 0),Vec(0.05, 0.50, 0.05))
objz = Rect(Vec(0, 0, 0),Vec(0.05, 0.05, 0.50))

setobject!(vis2[:world][:x], objx, matx)
setobject!(vis2[:world][:y], objy, maty)
setobject!(vis2[:world][:z], objz, matz)

setobject!(vis2[:body][:x], objx, matx)
setobject!(vis2[:body][:y], objy, maty)
setobject!(vis2[:body][:z], objz, matz)

setobject!(vis2[:vel][:x], objx, matx)
setobject!(vis2[:vel][:y], objy, maty)
setobject!(vis2[:vel][:z], objz, matz)

orientation = :MRP
rt = [0.04, 0.00, 0.00]
R = eval(orientation)(rt...)
@show R[:,3]
settransform!(vis2[:body], LinearMap(R))






@variables rt_[1:3]
rt_
R_ = Matrix(MRP(rt_...))
xb_ = R_[:,1]
yb_ = R_[:,2]
zb_ = R_[:,3]

Symbolics.jacobian(xb_, rt_)
Symbolics.jacobian(yb_, rt_)
Symbolics.jacobian(zb_, rt_)







# vis3 = Visualizer()
# open(vis3)
plot_surface!(vis3, model.env)
build_robot!(vis3, model, name=:test)
set_robot!(vis3, model, ref_traj.q[1], name=:test, offset=0.00)
θ = ref_traj.q[1][3]/π

q

Rotations.∇differential(R)

const ContactControl = Main

model = ContactControl.get_model("quadruped")
# model = ContactControl.get_model("hopper_2D")

# get trajectory
ref_traj = get_trajectory("quadruped", "gait1", load_type=:split_traj_alt, model=model)

# Sizes
nq = model.dim.q
nc = model.dim.c
nb = model.dim.b
nx = nq
ny = nb + 2nc
nz = ContactControl.num_var(model)
nθ = ContactControl.num_data(model)

# Indices
rz_ = ContactControl.RZLin(model, rand(nz,nz))
ix = rz_.ix
iy1 = rz_.iy1
iy2 = rz_.iy2
idyn = rz_.idyn
irst = rz_.irst
ibil = rz_.ibil

z0 = rand(nz)
θ0 = rand(nθ)
r0 = rand(nz)
rz0 = zeros(nz,nz)
rθ0 = zeros(nz,nθ)
model.res.rz!(rz0, z0, θ0)
model.res.rθ!(rθ0, z0, θ0)
vix = Vector(rz_.ix)
viy1 = Vector(rz_.iy1)
viy2 = Vector(rz_.iy2)
vidyn = Vector(rz_.idyn)
virst = Vector(rz_.irst)
vibil = Vector(rz_.ibil)
plot(Gray.(1e10*abs.(rz0[[vidyn; virst; vibil;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz0[[vidyn;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz0[[virst;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz0[[vibil;], [vix; viy1; viy2]])))

plot(Gray.(abs.(rz0[[vidyn; virst; vibil;], [vix; viy1; viy2]])))
plot(Gray.(abs.(rz0[[virst;], [viy1; ]])))



z = rand(nz)
θ = rand(nθ)
κ = 1e-4
κv = [κ]
rz1 = ContactControl.RZLin(model, rz0)
r1 = ContactControl.RLin(model, z0, θ0, r0, rz0, rθ0)

# Test rz!
ContactControl.rz!(rz1, z)
rz2 = rand(nz, nz)
model.linearized.rz!(rz2, z, rz0)

@test norm(rz1.Dx  - rz2[idyn, ix],  Inf) < 1e-10
@test norm(rz1.Dy1 - rz2[idyn, iy1], Inf) < 1e-10
@test norm(rz1.Rx  - rz2[irst, ix],  Inf) < 1e-10
@test norm(rz1.Ry1 - rz2[irst, iy1],  Inf) < 1e-10
@test norm(rz1.Ry2 - diag(rz2[irst, iy2]),  Inf) < 1e-10
@test norm(rz1.y2  - diag(rz2[ibil, iy1]),  Inf) < 1e-10
@test norm(rz1.y1  - diag(rz2[ibil, iy2]),  Inf) < 1e-10

# Test r!
r1 = ContactControl.RLin(model, z0, θ0, r0, rz0, rθ0)
ContactControl.r!(r1, z, θ, κ)
r2  = rand(nz)
model.linearized.r!(r2, z, θ, κv, z0, θ0, r0, rz0, rθ0)
@test norm(r2[r1.idyn] - r1.rdyn, Inf) < 1e-10
@test norm(r2[r1.irst] - r1.rrst, Inf) < 1e-10
@test norm(r2[r1.ibil] - r1.rbil, Inf) < 1e-10





using MeshCat
# using ModelingToolkit
using RigidBodyDynamics
using Symbolics
# import SymPy

# function setup(vis, n_dof, sympy=false)
#     robot = create_robot_kuka_iiwa_14(vis)
#     x = MechanismState{sympy ? SymPy.Sym : Num}(robot.mechanism)
#     q = configuration(x)
#     if sympy
#         for i in 1:n_dof
#             q[i] = SymPy.symbols("q_$i", real = true)
#         end
#     else
#         for i in 1:n_dof
#             q[i] = Num(Sym{Real}(Symbol("q$i")))
#         end
#     end
#     return x
# end

# vis = Visualizer()
# open(vis)
robot = create_robot_kuka_iiwa_14(vis) # This is a 7 DOF robot

n_dof = 2 # Number of degrees of freedom not == 0, set to something no more than 7


# using Symbolics
x = MechanismState{Num}(robot.mechanism)
q = configuration(x)
for i in 1:n_dof
    q[i] = Num(Sym{Real}(Symbol("q$i")))
end
@time M = mass_matrix(x); # Don't try to print it, it will take forever. Just printing the first element can take 20 seconds


julia> x = setup(vis, 3); @time M = mass_matrix(x);
  3.721925 seconds (7.73 M allocations: 462.482 MiB, 2.41% gc time, 99.20% compilation time)

julia> for n in 3:7; x = setup(vis, n); @time M = mass_matrix(x); end
  0.006976 seconds (79.61 k allocations: 3.028 MiB)
  0.007024 seconds (80.71 k allocations: 3.072 MiB)
  0.007367 seconds (81.56 k allocations: 3.106 MiB)
  0.007341 seconds (82.07 k allocations: 3.131 MiB)
  0.007315 seconds (82.66 k allocations: 3.157 MiB)

# In addition to `RigidBodyDynamics`, we'll be using the `StaticArrays`
# package, used throughout `RigidBodyDynamics`, which provides stack-allocated, fixed-size arrays:
using RigidBodyDynamics
using LinearAlgebra
using StaticArrays

# ## Creating a double pendulum `Mechanism`
# We're going to create a simple `Mechanism` that represents a
# [double pendulum](https://en.wikipedia.org/wiki/Double_pendulum).
# The `Mechanism` type can be thought of as an interconnection of rigid bodies and joints.
# We'll start by creating a 'root' rigid body, representing the
# fixed world, and using it to create a new `Mechanism`:

g = -9.81 # gravitational acceleration in z-direction
world = RigidBody{Float64}("world")
doublependulum = Mechanism(world; gravity = SVector(0, 0, g))


# Note that the `RigidBody` type is parameterized on the 'scalar type', here `Float64`.
# We'll now add a second body, called 'upper link', to the `Mechanism`.
# We'll attach it to the world with a revolute joint, with the $y$-axis as the
# axis of rotation. We'll start by creating a `SpatialInertia`, which stores the
# inertial properties of the new body:

axis = SVector(0., 1., 0.) # joint axis
I_1 = 0.333 # moment of inertia about joint axis
c_1 = -0.5 # center of mass location with respect to joint axis
m_1 = 1. # mass
frame1 = CartesianFrame3D("upper_link") # the reference frame in which the spatial inertia will be expressed
inertia1 = SpatialInertia(frame1,
    moment=I_1 * axis * axis',
    com=SVector(0, 0, c_1),
    mass=m_1)


# Note that the created `SpatialInertia` is annotated with the frame in which
# it is expressed (in the form of a `CartesianFrame3D`). This is a common theme
# among `RigidBodyDynamics` objects. Storing frame information with the data
# obviates the need for the complicated variable naming conventions that are
# used in some other libraries to disambiguate the frame in which quantities
# are expressed. It also enables automated reference frame checks.

# We'll now create the second body:
upperlink = RigidBody(inertia1)

# and a new revolute joint called 'shoulder':
shoulder = Joint("shoulder", Revolute(axis))

# Creating a `Joint` automatically constructs two new `CartesianFrame3D`
# objects: a frame directly before the joint, and one directly after. To attach
# the new body to the world by this joint, we'll have to specify where the
# frame before the joint is located on the parent body (here, the world):
before_shoulder_to_world = one(Transform3D,
    frame_before(shoulder), default_frame(world))

# Now we can attach the upper link to the world:
attach!(doublependulum, world, upperlink, shoulder,
    joint_pose = before_shoulder_to_world)


# which changes the tree representation of the `Mechanism`.
# We can attach the lower link in similar fashion:
l_1 = -1. # length of the upper link
I_2 = 0.333 # moment of inertia about joint axis
c_2 = -0.5 # center of mass location with respect to joint axis
m_2 = 1. # mass
inertia2 = SpatialInertia(CartesianFrame3D("lower_link"),
    moment=I_2 * axis * axis',
    com=SVector(0, 0, c_2),
    mass=m_2)
lowerlink = RigidBody(inertia2)
elbow = Joint("elbow", Revolute(axis))
before_elbow_to_after_shoulder = Transform3D(
    frame_before(elbow), frame_after(shoulder), SVector(0, 0, l_1))
attach!(doublependulum, upperlink, lowerlink, elbow,
    joint_pose = before_elbow_to_after_shoulder)


# Now our double pendulum `Mechanism` is complete.
# **Note**: instead of defining the `Mechanism` in this way, it is also
# possible to load in a [URDF](http://wiki.ros.org/urdf) file (an XML file
# format used in ROS), using the `parse_urdf` function, e.g.:
srcdir = dirname(pathof(RigidBodyDynamics))
urdf = joinpath(srcdir, "..", "test", "urdf", "Acrobot.urdf")
parse_urdf(urdf)

# ## The state of a `Mechanism`
# A `Mechanism` stores the joint/rigid body layout, but no state information.
# State information is separated out into a `MechanismState` object:
state = MechanismState(doublependulum)


# Let's first set the configurations and velocities of the joints:
set_configuration!(state, shoulder, 0.3)
set_configuration!(state, elbow, 0.4)
set_velocity!(state, shoulder, 1.)
set_velocity!(state, elbow, 2.);


# The joint configurations and velocities are stored as `Vector`s (denoted $q$
# and $v$ respectively in this package) inside the `MechanismState` object:
q = configuration(state)
v = velocity(state)


# ## Kinematics
# We are now ready to do kinematics. Here's how you transform a point at
# the origin of the frame after the elbow joint to world frame:
transform(state, Point3D(frame_after(elbow), zero(SVector{3})),
    default_frame(world))


# Other objects like `Wrench`es, `Twist`s, and `SpatialInertia`s can be transformed in similar fashion.
# You can also ask for the homogeneous transform to world:
transform_to_root(state, frame_after(elbow))

# Or a relative transform:
relative_transform(state, frame_after(elbow), frame_after(shoulder))

# and here's the center of mass of the double pendulum:
center_of_mass(state)

# ## Dynamics
# A `MechanismState` can also be used to compute quantities related to the
# dynamics of the `Mechanism`. Here we compute the mass matrix:
@btime mass_matrix(state)

# Note that there is also a zero-allocation version, `mass_matrix!`
# (the `!` at the end of a method is a Julia convention signifying that the
# function is 'in-place', i.e. modifies its input data).
# We can do inverse dynamics as follows (note again that there is a
# non-allocating version of this method as well):
v̇ = similar(velocity(state)) # the joint acceleration vector, i.e., the time
# derivative of the joint velocity vector v
v̇[shoulder][1] = 1
v̇[elbow][1] = 2
inverse_dynamics(state, v̇)


# ## Simulation
# Let's simulate the double pendulum for 5 seconds, starting from the state we
# defined earlier. For this, we can use the basic `simulate` function:
ts, qs, vs = simulate(state, 5., Δt = 1e-3);


# `simulate` returns a vector of times (`ts`) and associated joint
# configurations (`qs`) and velocities (`vs`). You can of course plot
# the trajectories using your favorite plotting package
# (see e.g. [Plots.jl](https://github.com/JuliaPlots/Plots.jl)).
# The [MeshCatMechanisms](https://github.com/JuliaRobotics/MeshCatMechanisms.jl)
# or [RigidBodyTreeInspector](https://github.com/rdeits/RigidBodyTreeInspector.jl)
# packages can also be used for 3D animation of the double pendulum in action.
# See also [RigidBodySim.jl](https://github.com/JuliaRobotics/RigidBodySim.jl)
# for a more full-fledge simulation environment.

# A lower level interface for simulation/ODE integration with more options is also available.
# Consult the documentation for more information.
# In addition, [RigidBodySim.jl](https://github.com/JuliaRobotics/RigidBodySim.jl)
# offers a more full-featured simulation environment.





srcdir = dirname(pathof(RigidBodyDynamics))
urdf = joinpath(srcdir, "..", "test", "urdf", "Acrobot.urdf")
doublependulum = parse_urdf(urdf)

state = MechanismState(doublependulum)

set_configuration!(state, shoulder, 0.3)
set_configuration!(state, elbow, 0.4)
set_velocity!(state, shoulder, 1.)
set_velocity!(state, elbow, 2.);

mass_matrix(state)


function code_gen(nq::Int)
    @variables q[1:nq]
    @variables w[1:nq]
    M = 2*Diagonal(q)
    M = Symbolics.simplify.(M)
    v = rand(Bool, nq) * q[1]
    v[1] += q[1]^2
    ∇v = Symbolics.jacobian(v, q, simplify=true)
    ∇Sv = Symbolics.sparsejacobian(v, q, simplify=true)
    ∇Sv_ = Symbolics.simplify(∇Sv.nzval)
    nz = length(∇Sv.nzval)
    w .*= 0.0
    w[1:nz, 1] = ∇Sv_
    a = 2.0 .* q

    expr = Dict{Symbol, Expr}()
    expr[:a1] = build_function(a, q)[1]
    expr[:a2] = build_function(a, q)[2]
    expr[:M1] = build_function(M, q)[1]
    expr[:M2] = build_function(M, q)[2]
    expr[:∇v1] = build_function(∇v, q)[1]
    expr[:∇v2] = build_function(∇v, q)[2]
    expr[:∇Sv1] = build_function(∇Sv_, q)[1]
    expr[:∇Sv2] = build_function(∇Sv_, q)[2]
    expr[:w1] = build_function(w, q)[1]
    expr[:w2] = build_function(w, q)[2]

    pattern = Dict{Symbol, AbstractSparseArray}()
    pattern[:∇Sv] = similar(∇Sv, Bool)
    return expr, pattern
end

T = Float64
nq = 40
expr, pattern = code_gen(nq)
a_fast1 = eval(expr[:a1])
a_fast2 = eval(expr[:a2])
mass_matrix_fast1 = eval(expr[:M1])
mass_matrix_fast2 = eval(expr[:M2])
jacobian_fast1 = eval(expr[:∇v1])
jacobian_fast2 = eval(expr[:∇v2])
sparsejacobian_fast1 = eval(expr[:∇Sv1])
sparsejacobian_fast2 = eval(expr[:∇Sv2])
vect_fast1 = eval(expr[:w1])
vect_fast2 = eval(expr[:w2])

expr[:a2]

q = zeros(T,nq)
qs = zeros(SVector{nq,T})
qsi = zeros(SizedVector{nq,T})

a = zeros(T,nq)
as = zeros(SVector{nq,T})
asi = zeros(SizedVector{nq,T})

M = zeros(T,nq,nq)
Ms = zeros(SMatrix{nq,nq,T,nq^2})
Msi = zeros(SizedMatrix{nq,nq,T})

∇v = zeros(T,nq,nq)
∇vs = zeros(SMatrix{nq,nq,T,nq^2})
∇vsi = zeros(SizedMatrix{nq,nq,T})

∇Sv = similar(pattern[:∇Sv], T)
nz = length(∇Sv.nzval)
∇Sv_ = zeros(T, nz)
∇Sv_s = zeros(SVector{nz,T})
∇Sv_si = zeros(SizedVector{nz,T})

w = zeros(T, nq)
ws = zeros(SVector{nq,T})
wsi = zeros(SizedVector{nq,T})




a_fast1(q)
a_fast1(qs)

a_fast2(a, q)
a_fast2(a, qs)
a_fast2(a, qsi)
a_fast2(asi, q)
a_fast2(asi, qs)
a_fast2(asi, qsi)

@benchmark a_fast1(q)
@benchmark a_fast1(qs)

@benchmark a_fast2(a, q)
@benchmark a_fast2(a, qs)
@benchmark a_fast2(a, qsi)
@benchmark a_fast2(asi, q)
@benchmark a_fast2(asi, qs)
@benchmark a_fast2(asi, qsi)


mass_matrix_fast1
mass_matrix_fast1(q)
mass_matrix_fast1(qs)
# mass_matrix_fast1(qsi)
mass_matrix_fast2(M, q)
mass_matrix_fast2(M, qs)
mass_matrix_fast2(M, qsi)
# mass_matrix_fast2(Ms, q)
# mass_matrix_fast2(Ms, qs)
# mass_matrix_fast2(Ms, qsi)
mass_matrix_fast2(Msi, q)
mass_matrix_fast2(Msi, qs)
mass_matrix_fast2(Msi, qsi)

@benchmark mass_matrix_fast1(q)
@benchmark mass_matrix_fast1(qs)
# @benchmark mass_matrix_fast1(qsi)
@benchmark mass_matrix_fast2(M, q)
@benchmark mass_matrix_fast2(M, qs)
@benchmark mass_matrix_fast2(M, qsi)
# @benchmark mass_matrix_fast2(Ms, q)
# @benchmark mass_matrix_fast2(Ms, qs)
# @benchmark mass_matrix_fast2(Ms, qsi)
@benchmark mass_matrix_fast2(Msi, q)
@benchmark mass_matrix_fast2(Msi, qs)
@benchmark mass_matrix_fast2(Msi, qsi)

jacobian_fast1(q)
jacobian_fast1(qs)

jacobian_fast2(∇v, q)
jacobian_fast2(∇v, qs)
jacobian_fast2(∇v, qsi)
jacobian_fast2(∇vsi, q)
jacobian_fast2(∇vsi, qs)
jacobian_fast2(∇vsi, qsi)

@benchmark jacobian_fast1(q)
@benchmark jacobian_fast1(qs)

@benchmark jacobian_fast2(∇v, q)
@benchmark jacobian_fast2(∇v, qs)
@benchmark jacobian_fast2(∇v, qsi)
@benchmark jacobian_fast2(∇vsi, q)
@benchmark jacobian_fast2(∇vsi, qs)
@benchmark jacobian_fast2(∇vsi, qsi)


sparsejacobian_fast1(q)
sparsejacobian_fast1(qs)

sparsejacobian_fast2(∇Sv_, q)
sparsejacobian_fast2(∇Sv_, qs)
sparsejacobian_fast2(∇Sv_, qsi)
sparsejacobian_fast2(∇Sv_si, q)
sparsejacobian_fast2(∇Sv_si, qs)
sparsejacobian_fast2(∇Sv_si, qsi)


@benchmark sparsejacobian_fast1(q)
@benchmark sparsejacobian_fast1(qs)

@benchmark sparsejacobian_fast2(∇Sv_, q)
@benchmark sparsejacobian_fast2(∇Sv_, qs)
@benchmark sparsejacobian_fast2(∇Sv_, qsi)
@benchmark sparsejacobian_fast2(∇Sv_si, q)
@benchmark sparsejacobian_fast2(∇Sv_si, qs)
@benchmark sparsejacobian_fast2(∇Sv_si, qsi)

vect_fast1(q)
vect_fast1(qs)

vect_fast2(w, q)
vect_fast2(w, qs)
vect_fast2(w, qsi)
vect_fast2(wsi, q)
vect_fast2(wsi, qs)
vect_fast2(wsi, qsi)

@benchmark vect_fast1(q)
@benchmark vect_fast1(qs)

@benchmark vect_fast2(w, q)
@benchmark vect_fast2(w, qs)
@benchmark vect_fast2(w, qsi)
@benchmark vect_fast2(wsi, q)
@benchmark vect_fast2(wsi, qs)
@benchmark vect_fast2(wsi, qsi)


# @variables q[1:nq]
# M = 2*Diagonal(q)
# M = Symbolics.simplify.(M)
# v = rand(Bool, nq) .* rand(nq) * q[1]
# v[1] += q[1]^2
# ∇v = Symbolics.jacobian(v, q, simplify=true)
# ∇Sv = Symbolics.sparsejacobian(v, q, simplify=true)
# ∇Sv_ = Symbolics.simplify(∇Sv.nzval)
#
# expr[:w2]



@variables x[1:10]
y = Diagonal(rand(Bool, 10) .* x.^2)

Symbolics.jacobian_sparsity(y, x)
