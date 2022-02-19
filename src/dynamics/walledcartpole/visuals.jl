function build_robot!(vis::Visualizer, model::WalledCartpole; name::Symbol=:WalledCartpole,
		r=0.05, α=1.0)
	r = convert(Float32, r)
	r_contact = convert(Float32, r * 1.5)

	wall_thickness = 0.15
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))
	wall_mat = MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0))

	default_background!(vis)

	link = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l), r/2)
	setobject!(vis[name][:robot]["link"], link, body_mat)
	vec_base = Vec(2r_contact, 4r, 1.5r)
	vec_shaft = Vec(3model.w + 2r_contact + 2wall_thickness, 0.5r, 0.5r)
	setobject!(vis[name][:robot]["linkbase"], Rect(-0.5 * vec_base, vec_base), body_mat)
	setobject!(vis[name][:robot]["shaft"], Rect(-0.5 * vec_shaft, vec_shaft), body_mat)
	setobject!(vis[name][:robot]["contact"], Sphere(Point3f0(0.0), r_contact), contact_mat)

	setobject!(vis[name][:env]["wall1"]["init"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness, wall_thickness, 0.75)), wall_mat)
	setobject!(vis[name][:env]["wall2"]["init"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness, wall_thickness, 0.75)), wall_mat)
	settransform!(vis[name][:env]["wall1"]["init"], Translation([-wall_thickness - model.w - 1.5 * r; -wall_thickness/2; 3r]))
	settransform!(vis[name][:env]["wall2"]["init"], Translation([model.w + 1.5 * r;                   -wall_thickness/2; 3r]))
	return nothing
end


function set_robot!(vis::Visualizer, model::WalledCartpole, q::AbstractVector;
		name::Symbol=:WalledCartpole)

	pt = [_kinematics(model, q, mode =  :tip)[1], 0.0, _kinematics(model, q, mode =  :tip)[2]]
    pb = [_kinematics(model, q, mode = :base)[1], 0.0, _kinematics(model, q, mode = :base)[2]]

	settransform!(vis[name][:robot]["link"], cable_transform(pt, pb))
	settransform!(vis[name][:robot]["linkbase"], Translation(pb))
	settransform!(vis[name][:robot]["contact"], Translation(pt))

	xw1 = q[3]
	xw2 = q[4]
	settransform!(vis[name][:env]["wall1"], Translation([xw1, 0.0, 0.0]))
	settransform!(vis[name][:env]["wall2"], Translation([xw2, 0.0, 0.0]))

	return nothing
end
