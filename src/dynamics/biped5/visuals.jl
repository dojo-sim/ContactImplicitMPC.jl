function plot_lines!(vis::Visualizer, model::Biped5, q::AbstractVector;
		r=0.035, offset=0.05, size=10, name::Symbol=:Biped5, col::Bool=true)
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

		k_toe_1 = kinematics_2(model, qi, body = :calf_1, mode = :ee)
		p_toe_1 = [k_toe_1[1], -offset, k_toe_1[2]] + p_shift

		k_toe_2 = kinematics_2(model, qi, body = :calf_2, mode = :ee)
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

function build_robot!(vis::Visualizer, model::Biped5; name::Symbol=:Biped5, r=0.035)
	r = convert(Float32, r)
	r_contact = convert(Float32, r*8/7)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0))

	default_background!(vis)

	torso = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_torso), r)
	setobject!(vis[name][:robot]["torso"], torso, body_mat)

	thigh_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_thigh1), r)
	setobject!(vis[name][:robot]["thigh1"], thigh_1, body_mat)

	calf_1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_calf1), r)
	setobject!(vis[name][:robot]["calf1"], calf_1, body_mat)

	thigh_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_thigh2), r)
	setobject!(vis[name][:robot]["thigh2"], thigh_2, body_mat)

	calf_2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l_calf2), r)
	setobject!(vis[name][:robot]["calf2"], calf_2, body_mat)

	setobject!(vis[name][:robot]["toe1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["toe2"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["knee1"], Sphere(Point3f0(0.0), r), body_mat)
	setobject!(vis[name][:robot]["knee2"], Sphere(Point3f0(0.0), r), body_mat)
	setobject!(vis[name][:robot]["hip"], Sphere(Point3f0(0.0), r), body_mat)
	setobject!(vis[name][:robot]["torso_top"], Sphere(Point3f0(0.0), r), body_mat)

	return nothing
end

function set_robot!(vis::Visualizer, model::Biped5, q::AbstractVector;
		name::Symbol=:Biped5, r=0.035)

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

	settransform!(vis[name][:robot]["thigh1"], cable_transform(p, p_thigh_1))
	settransform!(vis[name][:robot]["calf1"], cable_transform(p_thigh_1, p_calf_1))

	settransform!(vis[name][:robot]["thigh2"], cable_transform(p, p_thigh_2))
	settransform!(vis[name][:robot]["calf2"], cable_transform(p_thigh_2, p_calf_2))

	settransform!(vis[name][:robot]["torso"], cable_transform(p_torso,p))
	settransform!(vis[name][:robot]["toe1"], Translation(p_calf_1))
	settransform!(vis[name][:robot]["toe2"], Translation(p_calf_2))
	settransform!(vis[name][:robot]["knee1"], Translation(p_thigh_1))
	settransform!(vis[name][:robot]["knee2"], Translation(p_thigh_2))
	settransform!(vis[name][:robot]["hip"], Translation(p))
	settransform!(vis[name][:robot]["torso_top"], Translation(p_torso))

	return nothing
end

function contact_point(model::Biped5, q::AbstractVector)
	k_toe_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]]

	k_toe_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
	p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]]

	pc = [p_toe_1, p_toe_2]

	return pc
end
