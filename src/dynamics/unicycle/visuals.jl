function plot_lines!(vis::Visualizer, model::Unicycle, q::AbstractVector;
		offset=0.02, size=10, name::Symbol=:unicycle, col::Bool=true)
	p_shift = [0.0, 0.0, model.rw]

	# Point Traj
	body_point = Vector{Point{3,Float64}}()
	wheel_point = Vector{Point{3,Float64}}()

	for qi in q
		push!(body_point, Point(qi[1], -offset, qi[2])+p_shift)
		push!(wheel_point, Point(kinematics(model, qi)[1], -offset, kinematics(model, qi)[2])+p_shift)
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["/lines/body"], MeshCat.Line(body_point, orange_mat))
		setobject!(vis[name]["/lines/wheel"], MeshCat.Line(wheel_point, blue_mat))
	else
		setobject!(vis[name]["/lines/body"], MeshCat.Line(body_point, black_mat))
		setobject!(vis[name]["/lines/wheel"], MeshCat.Line(wheel_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Unicycle; name::Symbol=:Unicycle, r_sus=0.02, d=0.20, α=1.0)
	rw = convert(Float32, model.rw)
	r_sus = convert(Float32, r_sus)
	d = convert(Float32, d)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

    setobject!(vis[name][:robot][:body], HyperRectangle(-Vec(2d, d, 2d)/2,
		Vec(2d, d, 2d)), body_mat)
	setobject!(vis[name][:robot][:wheel][:tire], Cylinder(Point3f0(0), Point3f0(0, 0, d), rw), contact_mat)
	settransform!(vis[name][:robot][:wheel][:tire],
		compose(Translation([0, d/2, 0]),
			    LinearMap(RotX(π/2))))

	setobject!(vis[name][:robot][:wheel][:marker], Cylinder(Point3f0(0), Point3f0(0, 0, d/0.9), rw/5), body_mat)
	settransform!(vis[name][:robot][:wheel][:marker],
		compose(Translation([0, d/1.8, rw - rw/5]),
				LinearMap(RotX(π/2))))
    setobject!(vis[name][:robot][:suspension], Cylinder(Point3f0(0),
		Point3f0(0, 0, 1.0), r_sus), body_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::Unicycle, q::AbstractVector;
		name::Symbol=:Unicycle, r=0.10, r_sus=0.02, d=0.20)
	r = convert(Float32, r)
	r_sus = convert(Float32, r_sus)
	d = convert(Float32, d)
	p_shift = [0.0; 0.0; r]

    p_body = [q[1], 0.0, q[2]]
    p_wheel = [kinematics(model, q)[1], 0.0, kinematics(model, q)[2]]
	suspension_length = norm(p_body - p_wheel)

    settransform!(vis[name][:robot][:body], Translation(p_body))
    settransform!(vis[name][:robot][:wheel],
		compose(
			Translation(p_wheel),
			LinearMap(RotY(q[4]))),
			)
	settransform!(vis[name][:robot][:suspension],
		compose(
			cable_transform(p_wheel + 1e-10rand(3), p_body),
			LinearMap(Diagonal([1,1,suspension_length])))
		)
	return nothing
end

function contact_point(model::Unicycle, q::AbstractVector)
	p_wheel = [kinematics(model, q)[1], 0.0, kinematics(model, q)[2]]
	pc = [p_wheel]
	return pc
end
