module ContactImplicitMPC

using BenchmarkTools
using Colors
using CoordinateTransformations
using FFMPEG
using FileIO
using ForwardDiff
using GeometryBasics
using IfElse
using InteractiveUtils
using JLD2
using LinearAlgebra
using Logging
using QDLDL
using MeshCat
using MeshCatMechanisms
using Meshing
using MeshIO
using Parameters
using Plots
using Random
using Rotations
using SparseArrays
using StaticArrays
using SuiteSparse
using Symbolics
using Scratch
import Scratch: get_scratch!
using Test
using RoboDojo
import RoboDojo: LinearSolver, LUSolver, Model, ResidualMethods, Space, Disturbances, IndicesZ, InteriorPoint, EmptySolver, Policy, Trajectory, GradientTrajectory, InteriorPointOptions, IndicesOptimization, interior_point, interior_point_solve!, bilinear_violation, residual_violation, general_correction_term!, r!, rz!, rθ!, linear_solve!, lu_solver, empty_policy, empty_disturbances, friction_coefficients, SimulatorStatistics, SimulatorOptions, indices_θ, num_data, initialize_z!, initialize_θ!, indices_z, indices_θ, simulate!, policy, process!, Simulator
using DirectTrajectoryOptimization

# Utilities
include("utils.jl")

# Solver
include("solver/lu.jl") # sparse arrays
include("solver/ldl.jl")

include("solver/qr.jl")
include("solver/schur.jl")

# Environment
include("simulator/environment.jl")

# Dynamics
include("dynamics/model.jl")

# Simulator
include("simulation/index.jl")

# Simulator
include("simulation/contact_methods.jl")
include("simulation/simulation.jl")

include("dynamics/code_gen_dynamics.jl")
include("dynamics/fast_methods_dynamics.jl")

# Models
include("dynamics/quaternions.jl")
include("dynamics/mrp.jl")
include("dynamics/euler.jl")

include("dynamics/particle_2D/model.jl")
include("dynamics/particle/model.jl")
include("dynamics/hopper_2D/model.jl")
include("dynamics/hopper_3D/model.jl")
include("dynamics/quadruped/model.jl")
include("dynamics/flamingo/model.jl")
include("dynamics/pushbot/model.jl")
include("dynamics/walledcartpole/model.jl")
include("dynamics/centroidal_quadruped/model.jl")
include("dynamics/point_foot_quadruped/model.jl")

# Simulator
include("simulator/policy.jl")

include("simulator/disturbances.jl")
include("simulator/simulator.jl")

# Simulation
include("simulation/environments/flat.jl")
include("simulation/environments/piecewise.jl")
include("simulation/environments/quadratic.jl")
include("simulation/environments/slope.jl")
include("simulation/environments/sinusoidal.jl")
include("simulation/environments/stairs.jl")

include("simulation/residual_approx.jl")
include("simulation/code_gen_simulation.jl")

# Controller
include("controller/trajectory.jl")
include("controller/linearized_step.jl")
include("controller/implicit_dynamics.jl")
include("controller/objective.jl")
include("controller/linearized_solver.jl")
include("controller/newton.jl")
include("controller/newton_indices.jl")
include("controller/newton_residual.jl")
include("controller/newton_jacobian.jl")
include("controller/mpc_utils.jl")
include("controller/policy.jl")
include("controller/newton_structure_solver/methods.jl")

# Visuals
include("dynamics/visuals.jl")
include("dynamics/visual_utils.jl")
include("visuals.jl")

include("dynamics/particle_2D/visuals.jl")
include("dynamics/particle/visuals.jl")
include("dynamics/hopper_2D/visuals.jl")
include("dynamics/hopper_3D/visuals.jl")
include("dynamics/quadruped/visuals.jl")
include("dynamics/flamingo/visuals.jl")
include("dynamics/pushbot/visuals.jl")
include("dynamics/walledcartpole/visuals.jl")
include("dynamics/centroidal_quadruped/visuals.jl")
include("dynamics/point_foot_quadruped/visuals.jl")

export
    initial_conditions,
    World,
    LinearizedCone,
    NonlinearCone,
    Model,
    Dimensions,
    BaseMethods,
    DynamicsMethods,
    ResidualMethods,
    Environment,
    environment_2D,
    environment_3D,
    environment_2D_flat,
    environment_3D_flat,
    get_model,
    SparseStructure,
    LinearizedStep,
    bil_addition!,
    r_linearized!,
    rz_linearized!,
    ImplicitTrajectory,
    linearization!,
    implicit_dynamics!,
    TrackingObjective,
    TrackingVelocityObjective,
    second_order_cone_product,
    generate_base_expressions,
    RLin,
    RZLin,
    RθLin,
    ContactTraj,
    Simulation,
    num_var,
    num_data,
    get_simulation,
    get_trajectory,
    interior_point,
    InteriorPointOptions,
    interior_point_solve!,
    r!,
    rz!,
    rθ,
    generate_base_expressions,
    save_expressions,
    instantiate_base!,
    generate_dynamics_expressions,
    save_expressions,
    instantiate_dynamics!,
    environment_3D_flat,
    friction_dim,
    dim,
    sqrt_quat,
    cayley_map,
    L_multiply,
    R_multiply,
    R2,
    R3,
    FrictionCone,
    rotation,
    module_dir,
    open_loop_disturbances,
    disturbances,
    Disturbances,
    NoDisturbances,
    OpenLoopDisturbance,
    impulse_disturbances,
    ImpulseDisturbance,
    RandomDisturbance,
    open_loop_policy,
    policy,
    ci_mpc_policy,
    NewtonOptions,
    CIMPCOptions,
    SimulatorOptions,
    simulator,
    simulate!,
    generate_residual_expressions,
    instantiate_residual!,
    ϕ_func,
    tracking_error,
    repeat_ref_traj,
    Schur,
    schur_factorize!,
    schur_solve!,
    LinearSolver,
    LUSolver,
    lu_solver,
    ldl_solver,
    LDLSolver,
    factorize!,
    linear_solve!,
    IndicesOptimization,
    index_q2,
    index_γ1,
    index_b1,
    index_ψ1,
    index_s1,
    index_η1,
    index_s2,
    index_q0,
    index_q1,
    index_u1,
    index_w1,
    index_μ,
    index_h,
    index_dyn,
    index_imp,
    index_mdp,
    index_fri,
    index_bimp,
    index_bmdp,
    index_bfri,
    linearization_var_index,
    linearization_term_index,
    index_ort,
    index_soc,
    num_var,
    num_data,
    num_bilinear,
    index_equr,
    index_ortr,
    index_socr,
    visualize_meshrobot!,
    visualize_robot!,
    visualize_force!,
    visualize_disturbance!,
    visualize_payload!,
    process!,
    contact_trajectory,
    pack_z,
    pack_θ,
    generate_pusher_traj,
    update_friction_coefficient!
    quadratic_objective,
    particle,
    particle_2D,
    hopper_2D,
    hopper_3D,
    quadruped,
    quadruped_payload,
    flamingo,
    pushbot,
    walledcartpole,
    centroidal_quadruped,
    flat_3D_lc,
    flat_3D_nc,
    quadratic_bowl_3D_lc,
    flat_2D_lc,
    flat_2D_nc,
    slope1_2D_lc,
    sine2_2D_lc,
    stairs3_2D_lc,
    flat_3D_lc,
    sine2_3D_lc,
    sine1_2D_lc,
    piecewise1_2D_lc,
    sine3_2D_lc,
    slope_smooth_2D_lc,
    flat_2D_lc,
    flat_3D_lc,
    plot_surface!,
    stairs!,
    Simulator,
    set_trajectory!,
    initialize_z!

end # module
