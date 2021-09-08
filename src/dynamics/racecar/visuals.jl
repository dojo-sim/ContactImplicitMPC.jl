function plot_lines!(vis::Visualizer, model::Racecar12, q::AbstractVector;
		offset=0.02, size=10, name::Symbol=:racecar, col::Bool=true)

	# Point Traj
	body_point = Vector{Point{3,Float64}}()
	wheel_points = [Vector{Point{3,Float64}}() for i=1:4]

	for qi in q
		push!(body_point, Point(qi[1], qi[2], qi[3]))
		p1 = kinematics_2(model, qi, body = :wheel_hub1)
		p2 = kinematics_2(model, qi, body = :wheel_hub2)
		p3 = kinematics_2(model, qi, body = :wheel_hub3)
		p4 = kinematics_2(model, qi, body = :wheel_hub4)
		push!(wheel_points[1], Point(p1[1], p1[2], p1[3]))
		push!(wheel_points[2], Point(p2[1], p2[2], p2[3]))
		push!(wheel_points[3], Point(p3[1], p3[2], p3[3]))
		push!(wheel_points[4], Point(p4[1], p4[2], p4[3]))
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["lines/body"], MeshCat.Line(body_point, orange_mat))
		setobject!(vis[name]["lines/wheel1"], MeshCat.Line(wheel_points[1], blue_mat))
		setobject!(vis[name]["lines/wheel2"], MeshCat.Line(wheel_points[2], blue_mat))
		setobject!(vis[name]["lines/wheel3"], MeshCat.Line(wheel_points[3], blue_mat))
		setobject!(vis[name]["lines/wheel4"], MeshCat.Line(wheel_points[4], blue_mat))
	else
		setobject!(vis[name]["lines/body"], MeshCat.Line(body_point, black_mat))
		setobject!(vis[name]["lines/wheel1"], MeshCat.Line(wheel_points[1], black_mat))
		setobject!(vis[name]["lines/wheel2"], MeshCat.Line(wheel_points[2], black_mat))
		setobject!(vis[name]["lines/wheel3"], MeshCat.Line(wheel_points[3], black_mat))
		setobject!(vis[name]["lines/wheel4"], MeshCat.Line(wheel_points[4], black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Racecar12; name::Symbol=:Racecar12, rs=0.02, d=0.075, α=1.0)
	lb = convert(Float32, model.lb)
	wb = convert(Float32, model.wb)
	rw = convert(Float32, model.rw)
	rs = convert(Float32, rs)
	d = convert(Float32, d)
	orange_mat, blue_mat, black_mat = get_material(;α=1.0)


	default_background!(vis)

    setobject!(vis[name][:robot][:body], HyperRectangle(-Vec(lb+d, wb-3d, d)/2,
		Vec(lb+d, wb-3d, d)), black_mat)

	for i = 1:4
		setobject!(vis[name][:robot]["wheel$i"][:tire], Cylinder(Point3f0(0), Point3f0(0, 0, d), rw), orange_mat)
		settransform!(vis[name][:robot]["wheel$i"][:tire],
			compose(Translation([0, d/2, 0]),
				    LinearMap(RotX(π/2))))

		setobject!(vis[name][:robot]["wheel$i"][:marker], Cylinder(Point3f0(0), Point3f0(0, 0, d/0.9), rw/5), black_mat)
		settransform!(vis[name][:robot]["wheel$i"][:marker],
			compose(Translation([0, d/1.8, rw - rw/4]),
					LinearMap(RotX(π/2))))

		setobject!(vis[name][:robot]["suspension$i"], Cylinder(Point3f0(0),
			Point3f0(0, 0, 1.0), rs), blue_mat)
		@show i
	end
	return nothing
end

function set_robot!(vis::Visualizer, model::Racecar12, q::AbstractVector;
		name::Symbol=:Racecar12, r=0.10, rs=0.02, d = 0.075)
	r = convert(Float32, r)
	rs = convert(Float32, rs)
	ϵ = 1e-10rand(3)
	mirror = [2, 1, 4, 3]

	θ = q[3]
	lf = q[4]
	lr = q[5]

    p_body = [q[1], q[2], q[3]]
	p1 = kinematics_1(model, q, body = :wheel1)
	p2 = kinematics_1(model, q, body = :wheel2)
	p3 = kinematics_1(model, q, body = :wheel3)
	p4 = kinematics_1(model, q, body = :wheel4)
	hub1 = kinematics_2(model, q, body = :wheel_hub1)
	hub2 = kinematics_2(model, q, body = :wheel_hub2)
	hub3 = kinematics_2(model, q, body = :wheel_hub3)

    settransform!(vis[name][:robot][:body], compose(Translation(p_body), LinearMap(MRP(q[4:6]...))))

	for i = 1:4
		p = kinematics_1(model, q, body = Symbol("wheel$i"))
		p_mirror = kinematics_1(model, q, body = Symbol("wheel$(mirror[i])"))
		hub = kinematics_2(model, q, body = Symbol("wheel_hub$i"))
		α_rot = i <= 2 ? RotZ(q[7]) : RotZ(0.0)
		Δsus = p_mirror - p
		Δsus *= d / norm(Δsus)
		settransform!(vis[name][:robot]["wheel$i"],
			compose(
				Translation(hub),
				compose(
				compose(
				LinearMap(MRP(0.0, 0.0, q[6])),
				LinearMap(α_rot)),
				LinearMap(RotY(q[11+i])))
				))
		settransform!(vis[name][:robot]["suspension$i"],
			compose(
				cable_transform(hub + Δsus + ϵ, p + Δsus),
				LinearMap(Diagonal([1,1,q[7+i]])))
			)
	end
	return nothing
end

function contact_point(model::Racecar12, q::AbstractVector)
	p1 = kinematics_3(model, q, body = :wheel_contact1)
	p2 = kinematics_3(model, q, body = :wheel_contact2)
	p3 = kinematics_3(model, q, body = :wheel_contact3)
	p4 = kinematics_3(model, q, body = :wheel_contact4)
	pc = [p1, p2, p3, p4]
	return pc
end
