function plot_lines!(vis::Visualizer, model::Bicycle, q::AbstractVector;
		offset=0.02, size=10, name::Symbol=:bicycle, col::Bool=true)

	# Point Traj
	body_point = Vector{Point{3,Float64}}()
	wheel_points = [Vector{Point{3,Float64}}() for i=1:2]

	for qi in q
		push!(body_point, Point(qi[1], -offset, qi[2]))
		p_front = kinematics_2(model, qi, body = :front_hub)
		p_rear  = kinematics_2(model, qi, body = :rear_hub)
		push!(wheel_points[1], Point(p_front[1], -offset, p_front[2]))
		push!(wheel_points[2], Point(p_rear[1],  -offset, p_rear[2]))
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["/lines/body"], MeshCat.Line(body_point, orange_mat))
		setobject!(vis[name]["/lines/front_wheel"], MeshCat.Line(wheel_points[1], blue_mat))
		setobject!(vis[name]["/lines/rear_wheel"], MeshCat.Line(wheel_points[2], blue_mat))
	else
		setobject!(vis[name]["/lines/body"], MeshCat.Line(body_point, black_mat))
		setobject!(vis[name]["/lines/front_wheel"], MeshCat.Line(wheel_points[1], blue_mat))
		setobject!(vis[name]["/lines/rear_wheel"], MeshCat.Line(wheel_points[2], blue_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Bicycle; name::Symbol=:Bicycle, r_sus=0.02, d=0.20, α=1.0)
	lb = convert(Float32, model.lb)
	rw = convert(Float32, model.rw)
	r_sus = convert(Float32, r_sus)
	d = convert(Float32, d)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

    setobject!(vis[name][:robot][:body], HyperRectangle(-Vec(lb+d, d, d)/2,
		Vec(lb+d, d, d)), body_mat)

	setobject!(vis[name][:robot][:front_wheel][:tire], Cylinder(Point3f0(0), Point3f0(0, 0, d), rw), contact_mat)
	settransform!(vis[name][:robot][:front_wheel][:tire],
		compose(Translation([0, d/2, 0]),
			    LinearMap(RotX(π/2))))

	setobject!(vis[name][:robot][:front_wheel][:marker], Cylinder(Point3f0(0), Point3f0(0, 0, d/0.9), rw/5), body_mat)
	settransform!(vis[name][:robot][:front_wheel][:marker],
		compose(Translation([0, d/1.8, rw - rw/5]),
				LinearMap(RotX(π/2))))

	setobject!(vis[name][:robot][:rear_wheel][:tire], Cylinder(Point3f0(0), Point3f0(0, 0, d), rw), contact_mat)
	settransform!(vis[name][:robot][:rear_wheel][:tire],
		compose(Translation([0, d/2, 0]),
			    LinearMap(RotX(π/2))))

	setobject!(vis[name][:robot][:rear_wheel][:marker], Cylinder(Point3f0(0), Point3f0(0, 0, d/0.9), rw/5), body_mat)
	settransform!(vis[name][:robot][:rear_wheel][:marker],
		compose(Translation([0, d/1.8, rw - rw/5]),
				LinearMap(RotX(π/2))))

    setobject!(vis[name][:robot][:front_suspension], Cylinder(Point3f0(0),
		Point3f0(0, 0, 1.0), r_sus), body_mat)
    setobject!(vis[name][:robot][:rear_suspension], Cylinder(Point3f0(0),
		Point3f0(0, 0, 1.0), r_sus), body_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::Bicycle, q::AbstractVector;
		name::Symbol=:Bicycle, r=0.10, r_sus=0.02, d=0.20)
	r = convert(Float32, r)
	r_sus = convert(Float32, r_sus)
	d = convert(Float32, d)

	θ = q[3]
	lf = q[4]
	lr = q[5]

    p_body = [q[1], 0.0, q[2]]
	p_front = cast3d(kinematics_1(model, q, body = :front))
    p_rear = cast3d(kinematics_1(model, q, body = :rear))
	p_front_hub = cast3d(kinematics_2(model, q, body = :front_hub))
    p_rear_hub = cast3d(kinematics_2(model, q, body = :rear_hub))

    settransform!(vis[name][:robot][:body], compose(Translation(p_body), LinearMap(RotY(-θ))))
	settransform!(vis[name][:robot][:front_wheel],
		compose(
			Translation(p_front_hub),
			LinearMap(RotY(-q[6]))),
			)
	settransform!(vis[name][:robot][:rear_wheel],
		compose(
			Translation(p_rear_hub),
			LinearMap(RotY(-q[7]))),
			)
	settransform!(vis[name][:robot][:front_suspension],
		compose(
			cable_transform(p_front_hub + 1e-10rand(3), p_front),
			LinearMap(Diagonal([1,1,lf])))
		)
	settransform!(vis[name][:robot][:rear_suspension],
		compose(
			cable_transform(p_rear_hub + 1e-10rand(3), p_rear),
			LinearMap(Diagonal([1,1,lr])))
		)
	return nothing
end

function contact_point(model::Bicycle, q::AbstractVector)
	p_front_contact = cast3d(kinematics_3(model, q, body = :front_contact))
	p_rear_contact  = cast3d(kinematics_3(model, q, body = :rear_contact))
	pc = [p_front_contact, p_rear_contact]
	return pc
end
