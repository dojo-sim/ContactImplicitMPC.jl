include(joinpath(@__DIR__, "..", "dynamics", "biped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
# model = get_model("quadruped")
model_sim = get_model("biped", surf="sinusoidal")
model = get_model("biped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
# ref_traj = get_trajectory("biped", "gait5", load_type=:split_traj, model=model)
# ref_traj = get_trajectory("biped", "biped_gait (1)", load_type=:split_traj, model=model)
# ref_traj = get_trajectory("biped", "biped_gait (2)", load_type=:split_traj, model=model)
ref_traj = get_trajectory("biped", "biped_gait (3)", load_type=:split_traj, model=model)
# ref_traj = get_trajectory("biped", "biped_gait (5)", load_type=:split_traj, model=model)
visualize!(vis, model, ref_traj.q, Δt=h, name=:mpc)

ref_traj_copy = deepcopy(ref_traj)

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 800

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    # q = [Diagonal(1e-1 * [1.0, 0.01, 0.5, .15, .15, .15, .15, .01, .01]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [1.0, 0.01, 0.05, 1.5, 1.5, .15, .15, .0005, .0005]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [10; 1; 10; ones(nu-5); 10; 10]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        )
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

w_amp = [+0.02, -0.20]
sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    # d = random_disturbances(model, w_amp, H, h)
    # d = open_loop_disturbances([rand(model.dim.w) .* w_amp for i=1:H_sim]),
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

# @profiler status = simulate!(sim)
@time status = simulate!(sim)

l = 9
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i][l:l], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[l:l] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([b[1:nb] for b in sim.traj.b]*N_sample)...)', color=:red, linewidth=1.0)

visualize!(vis, model, sim.traj.q[1:N_sample:end], Δt=10*h/N_sample, name=:mpc)
visualize!(vis, model, sim.traj.q[1:N_sample:end], Δt=h, name=:mpc)
plot_lines!(vis, model, sim.traj.q[1:N_sample:end])
plot_surface!(vis, model_sim.env)


function plot_contact_force(vis::Visualizer, model::ContactDynamicsModel, traj::ContactTraj)
	nc = model.dim.c
	nu = model.dim.u
	nf = Int(nb / nc)
	ng = Int(nf/2)+1

	H = traj.H
	q = traj.q[3:end]
	γ = traj.γ
	b = traj.b

	pc = [[zeros(3)   for i=1:nc] for t=1:H]
	γc = [zeros(nc*ng)   for t=1:H]
	bc = [zeros(nc*ng)   for t=1:H]
	λc = [zeros(nc*ng)   for t=1:H]

	for t = 1:H
		k_toe_1 = kinematics_3(model, q[t], body = :foot_1, mode = :toe)
		p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]]

		k_heel_1 = kinematics_3(model, q[t], body = :foot_1, mode = :heel)
		p_heel_1 = [k_heel_1[1], 0.0, k_heel_1[2]]

		k_toe_2 = kinematics_3(model, q[t], body = :foot_2, mode = :toe)
		p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]]

		k_heel_2 = kinematics_3(model, q[t], body = :foot_2, mode = :heel)
		p_heel_2 = [k_heel_2[1], 0.0, k_heel_2[2]]

		pc[t] = [p_toe_1, p_heel_1, p_toe_2, p_heel_2]
		γc[t] = vcat([transpose(rotation(model.env, pc[t][i])) * [friction_mapping(model.env) * 0.0*b[t][(i-1) * nf .+ (1:nf)]; 1.0*γ[t][i]] for i = 1:nc]...) # TODO: make efficient
		bc[t] = vcat([transpose(rotation(model.env, pc[t][i])) * [friction_mapping(model.env) * 1.0*b[t][(i-1) * nf .+ (1:nf)]; 0.0*γ[t][i]] for i = 1:nc]...) # TODO: make efficient
		λc[t] = vcat([transpose(rotation(model.env, pc[t][i])) * [friction_mapping(model.env) * 1.0*b[t][(i-1) * nf .+ (1:nf)]; 1.0*γ[t][i]] for i = 1:nc]...) # TODO: make efficient

	end

	orange_mat, blue_mat, black_mat = get_line_material(10)
	anim = MeshCat.Animation(convert(Int, 100))

	point_x = [[Point(0.0,0.0,0.0), Point(1.0, 0.0, 0.0)] for i=1:nc]
	point_z = [[Point(0.0,0.0,0.0), Point(0.0, 0.0, 1.0)] for i=1:nc]
	for i = 1:nc
		setobject!(vis[:test][:impact]["$i"],   MeshCat.Line(point_z[i], blue_mat))
		setobject!(vis[:test][:friction]["$i"], MeshCat.Line(point_x[i], orange_mat))
	end

	p_shift = [0, -0.08, 0]
	for t = 1:H
		MeshCat.atframe(anim, t) do
			for i = 1:nc
				γi = γc[t][(i-1)*ng .+ (1:ng)]
				bi = bc[t][(i-1)*ng .+ (1:ng)]
				r =  rotation_3d(model.env, pc[t][i])

				trans = Translation(pc[t][i]+p_shift)
				scale_imp = LinearMap(r' * Diagonal([1.0, 1.0, norm(γi)]))
				transform = compose(trans, scale_imp)
				settransform!(vis[:test][:impact]["$i"], transform)

				trans = Translation(pc[t][i]+p_shift)
				scale_fri = LinearMap(r' * norm(bi))
				transform = compose(trans, scale_fri)
				settransform!(vis[:test][:friction]["$i"], transform)
			end
		end
	end
	MeshCat.setanimation!(vis, anim)
	return pc
end

get_contact_force(vis, model_sim, sim.traj)

a = 10
a = 10
a = 10


#
#
# function plot_impact(vis::Visualizer, model::ContactDynamicsModel, γ::AbstractVector;
# 		r=0.035, offset=0.05, size=10, name::Symbol=:biped, col::Bool=true)
# 	p_shift = [0.0, 0.0, r]
#
#     env = model.env
# 	orange_mat, blue_mat, black_mat = get_line_material(size)
#
# 	T = length(q)
# 	p_shift = [0.0; 0.0; r]
# 	for t = 1:T
# 		MeshCat.atframe(anim, t) do
# 			p = [q[t][1]; 0.0; q[t][2]] + p_shift
#
# 			k_torso = kinematics_1(model, q[t], body = :torso, mode = :ee)
# 			p_torso = [k_torso[1], 0.0, k_torso[2]] + p_shift
#
# 			k_thigh_1 = kinematics_1(model, q[t], body = :thigh_1, mode = :ee)
# 			p_thigh_1 = [k_thigh_1[1], 0.0, k_thigh_1[2]] + p_shift
#
# 			k_calf_1 = kinematics_2(model, q[t], body = :calf_1, mode = :ee)
# 			p_calf_1 = [k_calf_1[1], 0.0, k_calf_1[2]] + p_shift
#
# 			k_thigh_2 = kinematics_1(model, q[t], body = :thigh_2, mode = :ee)
# 			p_thigh_2 = [k_thigh_2[1], 0.0, k_thigh_2[2]] + p_shift
#
# 			k_calf_2 = kinematics_2(model, q[t], body = :calf_2, mode = :ee)
# 			p_calf_2 = [k_calf_2[1], 0.0, k_calf_2[2]] + p_shift
#
# 			k_toe_1 = kinematics_3(model, q[t], body = :foot_1, mode = :toe)
# 			p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]] + p_shift
#
# 			k_heel_1 = kinematics_3(model, q[t], body = :foot_1, mode = :heel)
# 			p_heel_1 = [k_heel_1[1], 0.0, k_heel_1[2]] + p_shift
#
# 			k_toe_2 = kinematics_3(model, q[t], body = :foot_2, mode = :toe)
# 			p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]] + p_shift
#
# 			k_heel_2 = kinematics_3(model, q[t], body = :foot_2, mode = :heel)
# 			p_heel_2 = [k_heel_2[1], 0.0, k_heel_2[2]] + p_shift
#
# 			settransform!(vis[name]["thigh1"], cable_transform(p, p_thigh_1))
# 			settransform!(vis[name]["calf1"], cable_transform(p_thigh_1, p_calf_1))
# 			settransform!(vis[name]["foot1"], cable_transform(p_toe_1, p_heel_1))
#
# 			settransform!(vis[name]["thigh2"], cable_transform(p, p_thigh_2))
# 			settransform!(vis[name]["calf2"], cable_transform(p_thigh_2, p_calf_2))
# 			settransform!(vis[name]["foot2"], cable_transform(p_toe_2, p_heel_2))
#
# 			settransform!(vis[name]["torso"], cable_transform(p_torso,p))
# 			settransform!(vis[name]["heel1"], Translation(p_heel_1))
# 			settransform!(vis[name]["heel2"], Translation(p_heel_2))
# 			settransform!(vis[name]["toe1"], Translation(p_toe_1))
# 			settransform!(vis[name]["toe2"], Translation(p_toe_2))
# 			settransform!(vis[name]["knee1"], Translation(p_thigh_1))
# 			settransform!(vis[name]["knee2"], Translation(p_thigh_2))
# 			settransform!(vis[name]["hip"], Translation(p))
# 			settransform!(vis[name]["torso_top"], Translation(p_torso))
# 		end
# 	end
#     return nothing
# end
#
#
#
# function plot_lines!(vis::Visualizer, model::Biped, q::AbstractVector;
# 		r=0.035, offset=0.05, size=10, name::Symbol=:biped, col::Bool=true)
# 	p_shift = [0.0, 0.0, r]
# 	orange_mat, blue_mat, black_mat = get_line_material(size)
#
# 	# Point Traj
# 	torso_point = Vector{Point{3,Float64}}()
# 	t1_point = Vector{Point{3,Float64}}()
# 	t2_point = Vector{Point{3,Float64}}()
# 	for qi in q
# 		k_torso = kinematics_1(model, qi, body = :torso, mode = :com)
# 		p_torso = [k_torso[1], -offset, k_torso[2]] + p_shift
#
# 		k_toe_1 = kinematics_3(model, qi, body = :foot_1, mode = :toe)
# 		p_toe_1 = [k_toe_1[1], -offset, k_toe_1[2]] + p_shift
#
# 		k_toe_2 = kinematics_3(model, qi, body = :foot_2, mode = :toe)
# 		p_toe_2 = [k_toe_2[1], -offset, k_toe_2[2]] + p_shift
#
# 		push!(torso_point, Point(p_torso...))
# 		push!(t1_point, Point(p_toe_1...))
# 		push!(t2_point, Point(p_toe_2...))
# 	end
#
# 	# Set lines
# 	if col
# 		setobject!(vis[name]["lines/torso"], MeshCat.Line(torso_point, orange_mat))
# 		setobject!(vis[name]["lines/toe1"], MeshCat.Line(t1_point, blue_mat))
# 		setobject!(vis[name]["lines/toe2"], MeshCat.Line(t2_point, blue_mat))
# 	else
# 		setobject!(vis[name]["lines/torso"], MeshCat.Line(torso_point, black_mat))
# 		setobject!(vis[name]["lines/toe1"], MeshCat.Line(t1_point, black_mat))
# 		setobject!(vis[name]["lines/toe2"], MeshCat.Line(t2_point, black_mat))
# 	end
# 	return nothing
# end
#
#
#
#
#
# function plot_friction()
#     return nothing
# end
#
# function plot_contact_force()
#     return nothing
# end
#
#
#
#
#
#
# a = 10
# a = 10
# a = 10
#
# # filename = "biped_sine"
# # MeshCat.convert_frames_to_video(
# #     "/home/simon/Downloads/$filename.tar",
# #     "/home/simon/Documents/$filename.mp4", overwrite=true)
# #
# # convert_video_to_gif(
# #     "/home/simon/Documents/$filename.mp4",
# #     "/home/simon/Documents/$filename.gif", overwrite=true)
#
# # const ContactControl = Main
#
# # qq = []
# # for q in ref_traj_copy.q
# #     for i = 1:N_sample
# #         push!(qq, q)
# #     end
# # end
# # plot(hcat(qq...)[1:model.dim.q, 1:end]',
# #     label = "", color = :black, width = 3.0)
# # plot!(hcat(sim.traj.q...)[1:model.dim.q, 1:100]',
# #     label = "", color = :cyan, width = 1.0, legend = :topleft)
# #
#
#
# # q_test = [ref_traj.q[1] + [0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.6, 0.6, ]]
# # visualize!(vis, model, q_test, Δt=10*h/N_sample, name=:mpc)
