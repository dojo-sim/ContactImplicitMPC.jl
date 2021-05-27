model = deepcopy(ContactControl.get_model("flamingo", surf = "flat"))
model.μ_world = 0.1
ref_traj = deepcopy(ContactControl.get_trajectory("flamingo", "gait0", load_type = :split_traj_alt))
update_friction_coefficient!(ref_traj, model)
H = ref_traj.H
h = ref_traj.h
for t = 1:H
	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
	@test norm(r) < 1.0e-4
end

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[2])

# simulator
sim = ContactControl.simulator(model, q0, q1, h, 40,
    p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u]),
    ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-5, solver = :lu_solver),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = ContactControl.simulate!(sim, verbose = false)
@test status
# @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.15
# end

include(joinpath(module_dir(), "src", "dynamics", "flamingo", "visuals.jl"))
vis = Visualizer()
render(vis)
anim = visualize_meshrobot!(vis, model, sim.traj)

qs = zeros(model.dim.q)

qs[4] = π / 10.0
qs[5] = -1.0 * qs[4]
qs[6] = qs[4]
qs[7] = -1.0 * qs[6]
qs[8] = 0.5 * π
qs[9] = 0.5 * π
k = kinematics(model, qs)[end]
qs[2] -= k
set_robot!(vis, model, qs, name = model_name(model))


# initial conditions
q0 = SVector{model.dim.q}(qs)
q1 = SVector{model.dim.q}(qs)

model.μ_world = 1.0
# simulator
sim = ContactControl.simulator(model, q0, q1, 0.01, 100,
    # p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u]),
    ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-5, solver = :lu_solver),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = ContactControl.simulate!(sim, verbose = false)
anim = visualize_robot!(vis, model, sim.traj)


# gravity compensating torque

idx_ineq = collect(1:0)
_num_var = model.dim.u + model.dim.c * dim(model.env)
_num_data = model.dim.q + 1

# residual
function _r!(r, z, θ, κ)
	u = z[1:model.dim.u]
	λ = z[model.dim.u .+ (1:model.dim.c * dim(model.env))]

	q = θ[1:model.dim.q]
	h = θ[model.dim.q .+ (1:1)]

	dynamics(model, h, q, q, u, zeros(model.dim.w), λ, q)
	nothing
end

@variables r_sym[1:model.dim.q]
@variables z_sym[1:_num_var]
@variables θ_sym[1:_num_data]
@variables κ_sym[1:1]

parallel = Symbolics.SerialForm()
_r!(r_sym, z_sym, θ_sym, κ_sym)
r_sym = simplify.(r_sym)
rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym,
	parallel = parallel)[2])
rz_exp = Symbolics.jacobian(r_sym, z_sym, simplify = true)
rθ_exp = Symbolics.jacobian(r_sym, θ_sym, simplify = true)
rz_sp = similar(rz_exp, Float64)
rθ_sp = similar(rθ_exp, Float64)
rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym,
	parallel = parallel)[2])
rθf! = eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym,
	parallel = parallel)[2])

rzf!(rz_sp, z, θ)

# options
opts = ContactControl.InteriorPointOptions(diff_sol = false)

# solver
z = 0.01 * randn(_num_var)
θ = [qs; h]
ip = ContactControl.interior_point(z, θ,
	idx_ineq = idx_ineq,
	r! = rf!, rz! = rzf!, rθ! = rθf!,
	rz = rz_sp,
	rθ = rθ_sp,
	opts = opts)

# solve
status = ContactControl.interior_point!(ip)
