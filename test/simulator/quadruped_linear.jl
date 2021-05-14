@testset "Simulator: Quadruped" begin
    # Reference trajectory
    model = deepcopy(ContactControl.get_model("quadruped", surf = "flat"))
    model.μ_world = 0.5

    ref_traj = deepcopy(ContactControl.get_trajectory("quadruped", "gait2", load_type = :split_traj_alt))
    ContactControl.update_friction_coefficient!(ref_traj, model)

    T = ref_traj.H
    h = ref_traj.h

    for t = 1:T
    	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
    	@test norm(r) < 1.0e-4
    end

    # initial conditions
    q0 = SVector{model.dim.q}(ref_traj.q[1])
    q1 = SVector{model.dim.q}(ref_traj.q[2])

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u]),
        ip_opts = ContactControl.InteriorPointOptions(
    		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-5, solver = :lu_solver),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    # simulate
    @test status = ContactControl.simulate!(sim, verbose = false)
    @show sim.traj.q[end][1:3]
    @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.025
end

# Reference trajectory
model = deepcopy(ContactControl.get_model("quadrupedlinear", surf = "flat"))
model.μ_world = 0.5
model.mb
model.mf
ref_traj = deepcopy(ContactControl.get_trajectory("quadrupedlinear", "gait1_mit_2.5percent", load_type = :split_traj_alt))
ContactControl.update_friction_coefficient!(ref_traj, model)

T = ref_traj.H
h = ref_traj.h

for t = 1:T
    r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
    @test norm(r) < 1.0e-5
end

plot(hcat(ref_traj.q...)[4:6, :]')

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 1
h_sim = h / N_sample
H_sim = H #4000 #3000

# barrier parameter
κ_mpc = 1.0e-4
H_mpc = 10

obj = TrackingVelocityObjective(H_mpc, model.dim,
    q = [Diagonal([[1.0; 0.1; 1.0]; [0.1; 0.1; 0.1]; 1.0e-3 * ones(12)]) for t = 1:H_mpc],
    u = [Diagonal(1.0e-3 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-8 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-8 * ones(model.dim.b)) for t = 1:H_mpc],
	v = [Diagonal(1.0e-5 * ones(model.dim.q)) ./ h^2 for t = 1:H_mpc])

p = linearized_mpc_policy(sim.traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        # solver = :ldl_solver,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions())

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
# mrp = MRP(RotX(0.1 * π))
# q1_ref[4:6] = [mrp.x; mrp.y; mrp.z]
# q0_ref[4:6] = [mrp.x; mrp.y; mrp.z]
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim_p = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
	# p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) .* [zeros(3); zeros(3); zeros(3); zeros(3)] for ut  in ref_traj.u]),
	# p = ContactControl.open_loop_policy([SVector{model.dim.u}(trim * h_sim) for ut in ref_traj.u]),
	p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6,
        diff_sol = false),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

time = @elapsed status = ContactControl.simulate!(sim_p)
# @elapsed status = ContactControl.simulate!(sim)
# @profiler status = ContactControl.simulate!(sim)
visualize!(vis, model, sim_p.traj.q, Δt = h_sim)
visualize!(vis, model, ref_traj.q, Δt = h_sim)

plot(hcat(ref_traj.q...)[1:3, :]', color = :black, width = 2.0)
plot!(hcat(sim.traj.q...)[1:3,:]', color = :red, width = 1.0)
plot(hcat(sim.traj.q...)[6 .+ (1:3),:]')


vis = Visualizer()
open(vis)

function visualize!(vis, model, q;
      r = 0.025, Δt = 0.1)

	default_background!(vis)

	setobject!(vis["torso"],
    	Rect(Vec(-model.l_torso, -model.w_torso, -0.05),Vec(2.0 * model.l_torso, 2.0 * model.w_torso, 0.05)),
    	MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	feet1 = setobject!(vis["feet1"], Sphere(Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet2 = setobject!(vis["feet2"], Sphere(Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet3 = setobject!(vis["feet3"], Sphere(Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet4 = setobject!(vis["feet4"], Sphere(Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	T = length(q)
	p_shift = [0.0, 0.0, r]

	for t = 1:T
		MeshCat.atframe(anim, t) do
			rot = MRP(q[t][4:6]...)

			p_torso = [q[t][1]; q[t][2]; q[t][3]] + p_shift
			p_foot1 = q[t][6 .+ (1:3)] + p_shift
			p_foot2 = q[t][9 .+ (1:3)] + p_shift
			p_foot3 = q[t][12 .+ (1:3)] + p_shift
			p_foot4 = q[t][15 .+ (1:3)] + p_shift

			settransform!(vis["torso"], compose(Translation(p_torso), LinearMap(MRP(q[t][4:6]...))))
			settransform!(vis["feet1"], Translation(p_foot1))
			settransform!(vis["feet2"], Translation(p_foot2))
			settransform!(vis["feet3"], Translation(p_foot3))
			settransform!(vis["feet4"], Translation(p_foot4))
		end
	end

	# settransform!(vis["/Cameras/default"],
	#     compose(Translation(0.0, 0.0, -1.0), LinearMap(RotZ(-pi / 2.0))))

	MeshCat.setanimation!(vis, anim)
end


plot_lines!(vis, model, sim.traj.q[1:25:end])
plot_surface!(vis, model.env, ylims=[0.3, -0.05])
anim = visualize_meshrobot!(vis, model, sim.traj, sample=5)

#
# model.res.rz!(model.spa.rz_sp, ref_traj.z[1], ref_traj.θ[1])
#
# show(sparse(model.spa.rz_sp))
