# function plot_lines!(vis::Visualizer, model::Hopper2D, q::AbstractVector;
# 		r_foot=0.04, offset=0.04, size=10, name::Symbol=:hopper_2D, col::Bool=true)
# 	p_shift = [0.0, 0.0, r_foot]
#
# 	kinematics(model::Hopper2D, q) = [q[1] + q[4] * sin(q[3]), q[2] - q[4] * cos(q[3])]
#
# 	# Point Traj
# 	top_point = Vector{Point{3,Float64}}()
# 	bot_point = Vector{Point{3,Float64}}()
#
# 	for qi in q
# 		push!(top_point, Point(qi[1], -offset, qi[2])+p_shift)
# 		push!(bot_point, Point(kinematics(model, qi)[1], -offset, kinematics(model, qi)[2])+p_shift)
# 	end
#
# 	# Set lines
# 	orange_mat, blue_mat, black_mat = get_line_material(size)
# 	if col
# 		setobject!(vis[name]["/lines/top"], MeshCat.Line(top_point, orange_mat))
# 		setobject!(vis[name]["/lines/bot"], MeshCat.Line(bot_point, blue_mat))
# 	else
# 		setobject!(vis[name]["/lines/top"], MeshCat.Line(top_point, black_mat))
# 		setobject!(vis[name]["/lines/bot"], MeshCat.Line(bot_point, black_mat))
# 	end
# 	return nothing
# end

function build_robot!(vis::Visualizer, model::InvertedPendulum; name::Symbol=:InvertedPendulum, r=0.05)
	n_leg = 100
	r = convert(Float32, r)
	r_contact = convert(Float32, r * 1.5)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0))

	default_background!(vis)

	link = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l), r)
	setobject!(vis[name][:robot]["link"], link, body_mat)

	setobject!(vis[name][:robot]["contact1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	# setobject!(vis[name][:robot]["contact2"], Sphere(Point3f0(0.0), r_contact), contact_mat)

	# settransform!(vis["/Cameras/default"],
	# 	compose(Translation(0.0, 0.5, -1.0), LinearMap(RotZ(-pi / 2.0))))
	return nothing
end

function add_walls!(vis::Visualizer, model::InvertedPendulum; name::Symbol=:InvertedPendulum, r=0.05)
	setobject!(vis[name][:env]["wall1"],
		Rect(Vec(0, 0, 0),Vec(1.0, 1.0, 1.5)),
		MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0)))

	setobject!(vis[name][:env]["wall2"],
		Rect(Vec(0, 0, 0),Vec(1.0, 1.0, 1.5)),
		MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0)))

	settransform!(vis[:InvertedPendulum][:env]["wall1"], Translation([0.25 + 1.5 * r; -0.5; 0.0]))
	settransform!(vis[:InvertedPendulum][:env]["wall2"], Translation([-1.25 - 1.5 * r; -0.5; 0.0]))
end

function set_robot!(vis::Visualizer, model::InvertedPendulum, q::AbstractVector;
		name::Symbol=:InvertedPendulum, r=0.04)
	r = convert(Float32, r)

	pee = [_kinematics(model, q, mode = :ee)[1], 0.0, _kinematics(model, q, mode = :ee)[2]]
    p1 = [_kinematics(model, q, mode = :d)[1], 0.0, _kinematics(model, q, mode = :d)[2]]

	settransform!(vis[name][:robot]["link"], cable_transform(zeros(3), pee))
	settransform!(vis[name][:robot]["contact1"], Translation(p1))
	# settransform!(vis[name][:robot]["contact2"], Translation(p2))

	return nothing
end
