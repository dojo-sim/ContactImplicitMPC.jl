function plot_lines!(vis::Visualizer, model::Flamingo, q::AbstractVector;
		r=0.005, offset=0.05, size=10, name::Symbol=:Flamingo, col::Bool=true)
	r_contact = r*8/7
	p_shift = [0.0, 0.0, r_contact]
	orange_mat, blue_mat, black_mat = get_line_material(size)

	# Point Traj
	torso_point = Vector{Point{3,Float64}}()
	t1_point = Vector{Point{3,Float64}}()
	t2_point = Vector{Point{3,Float64}}()
	for qi in q
		k_torso = kinematics_1(model, qi, body = :torso, mode = :com)
		p_torso = [k_torso[1], -offset, k_torso[2]] + p_shift

		k_toe_1 = kinematics_3(model, qi, body = :foot_1, mode = :toe)
		p_toe_1 = [k_toe_1[1], -offset, k_toe_1[2]] + p_shift

		k_toe_2 = kinematics_3(model, qi, body = :foot_2, mode = :toe)
		p_toe_2 = [k_toe_2[1], -offset, k_toe_2[2]] + p_shift

		push!(torso_point, Point(p_torso...))
		push!(t1_point, Point(p_toe_1...))
		push!(t2_point, Point(p_toe_2...))
	end

	# Set lines
	if col
		setobject!(vis[name][:lines][:torso], MeshCat.Line(torso_point, orange_mat))
		setobject!(vis[name][:lines][:toe1], MeshCat.Line(t1_point, blue_mat))
		setobject!(vis[name][:lines][:toe2], MeshCat.Line(t2_point, blue_mat))
	else
		setobject!(vis[name][:lines][:torso], MeshCat.Line(torso_point, black_mat))
		setobject!(vis[name][:lines][:toe1], MeshCat.Line(t1_point, black_mat))
		setobject!(vis[name][:lines][:toe2], MeshCat.Line(t2_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Flamingo; name::Symbol=:Flamingo, r=0.005, α=1.0)
	r = convert(Float32, r)
	r_contact = convert(Float32, r*8/7)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

	torso = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_torso), r)
	setobject!(vis[name][:robot]["torso"], torso, body_mat)

	thigh_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_thigh1), r)
	setobject!(vis[name][:robot]["thigh1"], thigh_1, body_mat)

	calf_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_calf1), r)
	setobject!(vis[name][:robot]["calf1"], calf_1, body_mat)

	foot_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_foot1 + model.d_foot1), r)
	setobject!(vis[name][:robot]["foot1"], foot_1, body_mat)

	thigh_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_thigh2), r)
	setobject!(vis[name][:robot]["thigh2"], thigh_2, body_mat)

	calf_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_calf2), r)
	setobject!(vis[name][:robot]["calf2"], calf_2, body_mat)

	foot_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_foot2 + model.d_foot2), r)
	setobject!(vis[name][:robot]["foot2"], foot_2, body_mat)

	setobject!(vis[name][:robot]["heel1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["heel2"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["toe1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["toe2"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["knee1"], Sphere(Point3f0(0.0), r), body_mat)
	setobject!(vis[name][:robot]["knee2"], Sphere(Point3f0(0.0), r), body_mat)
	setobject!(vis[name][:robot]["hip"], Sphere(Point3f0(0.0), r), body_mat)
	setobject!(vis[name][:robot]["torso_top"], Sphere(Point3f0(0.0), r), body_mat)

	return nothing
end

function set_robot!(vis::Visualizer, model::Flamingo, q::AbstractVector;
		name::Symbol=:Flamingo, r=0.005)

	r = convert(Float32, r)
	r_contact = convert(Float32, r*8/7)
	p_shift = [0.0; 0.0; r_contact]

	p = [q[1]; 0.0; q[2]] + p_shift

	k_torso = kinematics_1(model, q, body = :torso, mode = :ee)
	p_torso = [k_torso[1], 0.0, k_torso[2]] + p_shift

	k_thigh_1 = kinematics_1(model, q, body = :thigh_1, mode = :ee)
	p_thigh_1 = [k_thigh_1[1], 0.0, k_thigh_1[2]] + p_shift

	k_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	p_calf_1 = [k_calf_1[1], 0.0, k_calf_1[2]] + p_shift

	k_thigh_2 = kinematics_1(model, q, body = :thigh_2, mode = :ee)
	p_thigh_2 = [k_thigh_2[1], 0.0, k_thigh_2[2]] + p_shift

	k_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
	p_calf_2 = [k_calf_2[1], 0.0, k_calf_2[2]] + p_shift

	k_toe_1 = kinematics_3(model, q, body = :foot_1, mode = :toe)
	p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]] + p_shift

	k_heel_1 = kinematics_3(model, q, body = :foot_1, mode = :heel)
	p_heel_1 = [k_heel_1[1], 0.0, k_heel_1[2]] + p_shift

	k_toe_2 = kinematics_3(model, q, body = :foot_2, mode = :toe)
	p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]] + p_shift

	k_heel_2 = kinematics_3(model, q, body = :foot_2, mode = :heel)
	p_heel_2 = [k_heel_2[1], 0.0, k_heel_2[2]] + p_shift

	settransform!(vis[name][:robot]["thigh1"], cable_transform(p, p_thigh_1))
	settransform!(vis[name][:robot]["calf1"], cable_transform(p_thigh_1, p_calf_1))
	settransform!(vis[name][:robot]["foot1"], cable_transform(p_toe_1, p_heel_1))

	settransform!(vis[name][:robot]["thigh2"], cable_transform(p, p_thigh_2))
	settransform!(vis[name][:robot]["calf2"], cable_transform(p_thigh_2, p_calf_2))
	settransform!(vis[name][:robot]["foot2"], cable_transform(p_toe_2, p_heel_2))

	settransform!(vis[name][:robot]["torso"], cable_transform(p_torso,p))
	settransform!(vis[name][:robot]["heel1"], Translation(p_heel_1))
	settransform!(vis[name][:robot]["heel2"], Translation(p_heel_2))
	settransform!(vis[name][:robot]["toe1"], Translation(p_toe_1))
	settransform!(vis[name][:robot]["toe2"], Translation(p_toe_2))
	settransform!(vis[name][:robot]["knee1"], Translation(p_thigh_1))
	settransform!(vis[name][:robot]["knee2"], Translation(p_thigh_2))
	settransform!(vis[name][:robot]["hip"], Translation(p))
	settransform!(vis[name][:robot]["torso_top"], Translation(p_torso))

	return nothing
end

function contact_point(model::Flamingo, q::AbstractVector)
	k_toe_1 = kinematics_3(model, q, body = :foot_1, mode = :toe)
	p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]]

	k_heel_1 = kinematics_3(model, q, body = :foot_1, mode = :heel)
	p_heel_1 = [k_heel_1[1], 0.0, k_heel_1[2]]

	k_toe_2 = kinematics_3(model, q, body = :foot_2, mode = :toe)
	p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]]

	k_heel_2 = kinematics_3(model, q, body = :foot_2, mode = :heel)
	p_heel_2 = [k_heel_2[1], 0.0, k_heel_2[2]]

	pc = [p_toe_1, p_heel_1, p_toe_2, p_heel_2]
	return pc
end

function build_meshrobot!(vis::Visualizer, model::Flamingo; name::Symbol=:Flamingo, α=1.0)
	default_background!(vis)

	# urdf = joinpath(@__DIR__, "mesh", "flamingo.urdf")
	# mechanism = MeshCatMechanisms.parse_urdf(urdf)
	# mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf, package_path=[@__DIR__]), vis[name])
	urdf = joinpath(@__DIR__, "mesh", "flamingo.urdf")
	package_path = @__DIR__
	build_meshrobot!(vis, model, urdf, package_path; name=name, α=α)
end

function convert_config(model::Flamingo, q::AbstractVector)
    _q = zeros(7)
    _q[1] = q[3]
    _q[2] = q[4] - _q[1]
    _q[4] = (q[5] - _q[2] - _q[1])
    _q[6] = q[8] - _q[4] - _q[2] - _q[1] - 0.5 * π

    _q[3] = q[6] - _q[1]
    _q[5] = q[7] - _q[3] - _q[1]
    _q[7] = q[9] - _q[5] - _q[3] - _q[1] - 0.5 * π

    return _q, compose(Translation(q[1], 0.0, q[2] + 0.02), LinearMap(RotZ(π)))
end


function flamingo_ghost!(vis, sim, surf; lines = true)
	x_mean = 1.5 #0.5 * (sim.traj.q[1][1] + sim.traj.q[end][1])
	shift_traj = deepcopy(sim.traj)
	sim.traj.q[end][1]
	for t = 1:length(shift_traj.q)
		shift_traj.q[t][1] -= x_mean
	end

	# surf_shift = x -> x[3] - sim.model.env.surf(x[1] + x_mean)[1]
	# plot_surface!(vis, surf_shift, xlims=[-4.5, 4.5], ylims = [-0.5, 0.5],
	# 	col = (0.9, 0.9, 0.9), α = 1.0)
	x_range = range(-3.0, stop = 3.0, length = 100)
	surface_points = [Point(x, 0.0, surf(x + x_mean)) for x in x_range]
	surface_mat = LineBasicMaterial(color=RGBA(0.0, 0.0, 0.0, 1.0), linewidth=5)
	setobject!(vis[:lines][:surface], MeshCat.Line(surface_points, surface_mat))

	lines && plot_lines!(vis, model, shift_traj.q[1:N_sample:end], offset=-0.25, size = 5)

	settransform!(vis["/Cameras/default"],
			compose(Translation(0.0, -95.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 50)

	t = 1
	mvis1 = build_meshrobot!(vis, sim.model; name=:Flamingo1, α = 0.20)
	set_meshrobot!(vis, mvis1, sim.model, shift_traj.q[t], name = :Flamingo1)

	t = 1175
	mvis2 = build_meshrobot!(vis, sim.model; name=:Flamingo2, α = 0.40)
	set_meshrobot!(vis, mvis2, sim.model, shift_traj.q[t], name = :Flamingo2)

	t = 2565
	mvis3 = build_meshrobot!(vis, sim.model; name=:Flamingo3, α = 0.60)
	set_meshrobot!(vis, mvis3, sim.model, shift_traj.q[t], name = :Flamingo3)

	t = 3750
	mvis4 = build_meshrobot!(vis, sim.model; name=:Flamingo4, α = 0.80)
	set_meshrobot!(vis, mvis4, sim.model, shift_traj.q[t], name = :Flamingo4)

	t = shift_traj.H
	mvis5 = build_meshrobot!(vis, sim.model; name=:Flamingo5, α = 1.0)
	set_meshrobot!(vis, mvis5, sim.model, shift_traj.q[t], name = :Flamingo5)
end
