function plot_lines!(vis::Visualizer, model::Quadruped, q::AbstractVector;
		r=0.0205, offset=0.05, size=10, name::Symbol=:Quadruped, col::Bool=true)
	p_shift = [0.0, 0.0, r]
	orange_mat, blue_mat, black_mat = get_line_material(size)

	# Point Traj
	torso_point = Vector{Point{3,Float64}}()
	f1_point = Vector{Point{3,Float64}}()
	f2_point = Vector{Point{3,Float64}}()
	f3_point = Vector{Point{3,Float64}}()
	f4_point = Vector{Point{3,Float64}}()

	for qi in q
		k_torso = kinematics_1(model, qi, body = :torso, mode = :com)
		p_torso = [k_torso[1], -offset, k_torso[2]] + p_shift

		k_calf_1 = kinematics_2(model, qi, body = :calf_1, mode = :ee)
		p_calf_1 = [k_calf_1[1], -offset, k_calf_1[2]] + p_shift

		k_calf_2 = kinematics_2(model, qi, body = :calf_2, mode = :ee)
		p_calf_2 = [k_calf_2[1], -offset, k_calf_2[2]] + p_shift

		k_calf_3 = kinematics_3(model, qi, body = :calf_3, mode = :ee)
		p_calf_3 = [k_calf_3[1], -offset, k_calf_3[2]] + p_shift

		k_calf_4 = kinematics_3(model, qi, body = :calf_4, mode = :ee)
		p_calf_4 = [k_calf_4[1], -offset, k_calf_4[2]] + p_shift

		push!(torso_point, Point(p_torso...))
		push!(f1_point, Point(p_calf_1...))
		push!(f2_point, Point(p_calf_2...))
		push!(f3_point, Point(p_calf_3...))
		push!(f4_point, Point(p_calf_4...))
	end

	# Set lines
	if col
		setobject!(vis[name]["lines/torso"], MeshCat.Line(torso_point, orange_mat))
		setobject!(vis[name]["lines/foot1"], MeshCat.Line(f1_point, blue_mat))
		setobject!(vis[name]["lines/foot2"], MeshCat.Line(f2_point, blue_mat))
		setobject!(vis[name]["lines/foot3"], MeshCat.Line(f3_point, blue_mat))
		setobject!(vis[name]["lines/foot4"], MeshCat.Line(f4_point, blue_mat))
	else
		setobject!(vis[name]["lines/torso"], MeshCat.Line(torso_point, black_mat))
		setobject!(vis[name]["lines/foot1"], MeshCat.Line(f1_point, black_mat))
		setobject!(vis[name]["lines/foot2"], MeshCat.Line(f2_point, black_mat))
		setobject!(vis[name]["lines/foot3"], MeshCat.Line(f3_point, black_mat))
		setobject!(vis[name]["lines/foot4"], MeshCat.Line(f4_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Quadruped; name::Symbol=:Quadruped, r=0.0205)
	r = convert(Float32, r)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0))

	default_background!(vis)

	torso = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_torso),
		convert(Float32, 0.035))
	setobject!(vis[name][:robot]["torso"], torso, body_mat)

	thigh_1 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_thigh1),
		convert(Float32, 0.0175))
	setobject!(vis[name][:robot]["thigh1"], thigh_1, body_mat)

	calf_1 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_calf1),
		convert(Float32, 0.0125))
	setobject!(vis[name][:robot]["leg1"], calf_1, body_mat)

	thigh_2 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_thigh2),
		convert(Float32, 0.0175))
	setobject!(vis[name][:robot]["thigh2"], thigh_2, body_mat)

	calf_2 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_calf2),
		convert(Float32, 0.0125))
	setobject!(vis[name][:robot]["leg2"], calf_2, body_mat)

	thigh_3 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_thigh3),
		convert(Float32, 0.0175))
	setobject!(vis[name][:robot]["thigh3"], thigh_3, body_mat)

	calf_3 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_calf3),
		convert(Float32, 0.0125))
	setobject!(vis[name][:robot]["leg3"], calf_3, body_mat)

	thigh_4 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_thigh4),
		convert(Float32, 0.0175))
	setobject!(vis[name][:robot]["thigh4"], thigh_4, body_mat)

	calf_4 = GeometryBasics.Cylinder(Point3f0(0),
		Point3f0(0.0, 0.0, model.l_calf4),
		convert(Float32, 0.0125))
	setobject!(vis[name][:robot]["leg4"], calf_4, body_mat)

	hip1 = setobject!(vis[name][:robot]["hip1"],
		GeometryBasics.Sphere(Point3f0(0), convert(Float32, 0.035)), body_mat)
	hip2 = setobject!(vis[name][:robot]["hip2"],
		GeometryBasics.Sphere(Point3f0(0), convert(Float32, 0.035)), body_mat)
	knee1 = setobject!(vis[name][:robot]["knee1"],
		GeometryBasics.Sphere(Point3f0(0), r), body_mat)
	knee2 = setobject!(vis[name][:robot]["knee2"],
		GeometryBasics.Sphere(Point3f0(0), r), body_mat)
	knee3 = setobject!(vis[name][:robot]["knee3"],
		GeometryBasics.Sphere(Point3f0(0), r), body_mat)
	knee4 = setobject!(vis[name][:robot]["knee4"],
		GeometryBasics.Sphere(Point3f0(0), r), body_mat)

	feet1 = setobject!(vis[name][:robot]["feet1"],
		GeometryBasics.Sphere(Point3f0(0), r), contact_mat)
	feet2 = setobject!(vis[name][:robot]["feet2"],
		GeometryBasics.Sphere(Point3f0(0), r), contact_mat)
	feet3 = setobject!(vis[name][:robot]["feet3"],
		GeometryBasics.Sphere(Point3f0(0), r), contact_mat)
	feet4 = setobject!(vis[name][:robot]["feet4"],
		GeometryBasics.Sphere(Point3f0(0), r), contact_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::Quadruped, q::AbstractVector;
		name::Symbol=:Quadruped, r=0.0205, offset=0.00)

	r = convert(Float32, r)
	p_shift = [0.0; offset; r]

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


	k_thigh_3 = kinematics_2(model, q, body = :thigh_3, mode = :ee)
	p_thigh_3 = [k_thigh_3[1], 0.0, k_thigh_3[2]] + p_shift

	k_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :ee)
	p_calf_3 = [k_calf_3[1], 0.0, k_calf_3[2]] + p_shift

	k_thigh_4 = kinematics_2(model, q, body = :thigh_4, mode = :ee)
	p_thigh_4 = [k_thigh_4[1], 0.0, k_thigh_4[2]] + p_shift

	k_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :ee)
	p_calf_4 = [k_calf_4[1], 0.0, k_calf_4[2]] + p_shift

	settransform!(vis[name][:robot]["thigh1"], cable_transform(p, p_thigh_1))
	settransform!(vis[name][:robot]["leg1"], cable_transform(p_thigh_1, p_calf_1))
	settransform!(vis[name][:robot]["thigh2"], cable_transform(p, p_thigh_2))
	settransform!(vis[name][:robot]["leg2"], cable_transform(p_thigh_2, p_calf_2))
	settransform!(vis[name][:robot]["thigh3"], cable_transform(p_torso, p_thigh_3))
	settransform!(vis[name][:robot]["leg3"], cable_transform(p_thigh_3, p_calf_3))
	settransform!(vis[name][:robot]["thigh4"], cable_transform(p_torso, p_thigh_4))
	settransform!(vis[name][:robot]["leg4"], cable_transform(p_thigh_4, p_calf_4))
	settransform!(vis[name][:robot]["torso"], cable_transform(p, p_torso))
	settransform!(vis[name][:robot]["hip1"], MeshCat.Translation(p))
	settransform!(vis[name][:robot]["hip2"], MeshCat.Translation(p_torso))
	settransform!(vis[name][:robot]["knee1"], MeshCat.Translation(p_thigh_1))
	settransform!(vis[name][:robot]["knee2"], MeshCat.Translation(p_thigh_2))
	settransform!(vis[name][:robot]["knee3"], MeshCat.Translation(p_thigh_3))
	settransform!(vis[name][:robot]["knee4"], MeshCat.Translation(p_thigh_4))
	settransform!(vis[name][:robot]["feet1"], MeshCat.Translation(p_calf_1))
	settransform!(vis[name][:robot]["feet2"], MeshCat.Translation(p_calf_2))
	settransform!(vis[name][:robot]["feet3"], MeshCat.Translation(p_calf_3))
	settransform!(vis[name][:robot]["feet4"], MeshCat.Translation(p_calf_4))

	return nothing
end

function contact_point(model::Quadruped, q::AbstractVector)
	k_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	p_calf_1 = [k_calf_1[1], 0.0, k_calf_1[2]]

	k_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
	p_calf_2 = [k_calf_2[1], 0.0, k_calf_2[2]]

	k_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :ee)
	p_calf_3 = [k_calf_3[1], 0.0, k_calf_3[2]]

	k_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :ee)
	p_calf_4 = [k_calf_4[1], 0.0, k_calf_4[2]]

	pc = [p_calf_1, p_calf_2, p_calf_3, p_calf_4]
	return pc
end

function build_meshrobot!(vis::Visualizer, model::Quadruped; name::Symbol=:Quadruped, α=1.0)
	default_background!(vis)

	urdf = joinpath(@__DIR__, "mesh", "a1.urdf")
	package_path = @__DIR__
	build_meshrobot!(vis, model, urdf, package_path; name=name, α=α)
end



function convert_config(model::Quadruped, q::AbstractVector)
    # Quadruped configuration
    # 1.  position long axis +x
    # 2.  position long axis +z
    # 3.  trunk rotation along -y
    # 4.  back  left  shoulder rotation along -y
    # 5.  back  left  elbow    rotation along -y
    # 6.  back  right shoulder rotation along -y
    # 7.  back  right elbow    rotation along -y
    # 8.  front right shoulder rotation along -y
    # 9.  front right elbow    rotation along -y
	# 10. front left  shoulder rotation along -y
	# 11. front left  elbow    rotation along -y

    # URDF configuration
    # 1.  front right clavicle rotation along +x
    # 2.  front left  clavicle rotation along +x
    # 3.  back  right clavicle rotation along +x
    # 4.  back  left  clavicle rotation along +x
    # 5.  front right shoulder rotation along +y
    # 6.  front left  shoulder rotation along +y
    # 7.  back  right shoulder rotation along +y
    # 8.  back  left  shoulder rotation along +y
    # 9.  front right elbow    rotation relative to shoulder along +y
    # 10. front left  elbow    rotation relative to shoulder along +y
    # 11. back  right elbow    rotation relative to shoulder along +y
    # 12. back  left  elbow    rotation relative to shoulder along +y
    q_ = zeros(12)
    x, z, θ = q[1:3]
    q_[5]  = -q[10] + θ - π/2
    q_[6]  = -q[8] + θ - π/2
	q_[7]  = -q[4] + θ - π/2
    q_[8]  = -q[6] + θ - π/2
    q_[9]  = +q[10]-q[11]
    q_[10] = +q[8]-q[9]
	q_[11] = +q[4]-q[5]
    q_[12] = +q[6]-q[7]

	r_foot = 0.02
	trunk_length = 2*0.183
	trunk_width = 2*0.132
    tform = compose(
		Translation(x , trunk_width/2, z+r_foot),
		compose(
			LinearMap(AngleAxis(π/2-θ, 0, 1.0, 0)),
			Translation(trunk_length/2, 0,0)
			)
		)
    return q_, tform
end
