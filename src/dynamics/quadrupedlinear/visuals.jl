function plot_lines!(vis::Visualizer, model::QuadrupedLinear12, q::AbstractVector;
		r=0.0205, offset=0.05, size=10, name::Symbol=:QuadrupedLinear12, col::Bool=true, α::Real=1.0)
	p_shift = [0.0, 0.0, r]
	orange_mat, blue_mat, black_mat = get_line_material(size, α=α)
	black_mat.color = RGBA(0.0,153/256,153/256,α)
	black_mat.color = RGBA(209/256,0/256,209/256,α)
	black_mat.color = RGBA(136/256,008/256,008/256,α)
	black_mat.color = RGBA(170/256,074/256,068/256,α)

	# Point Traj
	torso_point = Vector{Point{3,Float64}}()
	f1_point = Vector{Point{3,Float64}}()
	f2_point = Vector{Point{3,Float64}}()
	f3_point = Vector{Point{3,Float64}}()
	f4_point = Vector{Point{3,Float64}}()

	for qi in q
		p_torso = q[1:3] + p_shift
		p_foot_1 = q[6  .+ [1:3]] + p_shift
		p_foot_2 = q[9  .+ [1:3]] + p_shift
		p_foot_3 = q[12 .+ [1:3]] + p_shift
		p_foot_4 = q[15 .+ [1:3]] + p_shift
		push!(torso_point, Point(p_torso...))
		push!(f1_point, Point(p_foot_1...))
		push!(f2_point, Point(p_foot_2...))
		push!(f3_point, Point(p_foot_3...))
		push!(f4_point, Point(p_foot_4...))
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


function build_robot!(vis::Visualizer, model::QuadrupedLinear12; name::Symbol=:QuadrupedLinear12, r=0.0205, α=1.0)
	r = convert(Float32, r)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

	setobject!(vis[name][:robot][:torso],
    	Rect(Vec(-model.l_torso, -model.w_torso, -0.05),
			Vec(2.0 * model.l_torso, 2.0 * model.w_torso, 0.05)),
		body_mat)

	foot1 = setobject!(vis[name][:robot][:foot_1], Sphere(Point3f0(0), r), contact_mat)
	foot2 = setobject!(vis[name][:robot][:foot_2], Sphere(Point3f0(0), r), contact_mat)
	foot3 = setobject!(vis[name][:robot][:foot_3], Sphere(Point3f0(0), r), contact_mat)
	foot4 = setobject!(vis[name][:robot][:foot_4], Sphere(Point3f0(0), r), contact_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::QuadrupedLinear12, q::AbstractVector;
		name::Symbol=:QuadrupedLinear12, r=0.0205, offset=0.00)

	r = convert(Float32, r)
	p_shift = [0.0; offset; r]

	rot = MRP(q[4:6]...)

	p_torso = [q[1]; q[2]; q[3]] + p_shift
	p_foot1 = q[6 .+ (1:3)] + p_shift
	p_foot2 = q[9 .+ (1:3)] + p_shift
	p_foot3 = q[12 .+ (1:3)] + p_shift
	p_foot4 = q[15 .+ (1:3)] + p_shift

	settransform!(vis[name][:robot][:torso], compose(Translation(p_torso), LinearMap(rot)))
	settransform!(vis[name][:robot][:foot_1], Translation(p_foot1))
	settransform!(vis[name][:robot][:foot_2], Translation(p_foot2))
	settransform!(vis[name][:robot][:foot_3], Translation(p_foot3))
	settransform!(vis[name][:robot][:foot_4], Translation(p_foot4))
	return nothing
end

function contact_point(model::QuadrupedLinear12, q::AbstractVector)
	p_foot_1 = q[6  .+ [1:3]]
	p_foot_2 = q[9  .+ [1:3]]
	p_foot_3 = q[12 .+ [1:3]]
	p_foot_4 = q[15 .+ [1:3]]

	pc = [p_foot_1, p_foot_2, p_foot_3, p_foot_4]
	return pc
end

function build_meshrobot!(vis::Visualizer, model::QuadrupedLinear12; name::Symbol=:QuadrupedLinear12, α=1.0)
	default_background!(vis)

	urdf = joinpath(@__DIR__, "mesh", "a1.urdf")
	package_path = @__DIR__
	build_meshrobot!(vis, model, urdf, package_path; name=name, α=α)
end

function convert_config(model::QuadrupedLinear12, q::AbstractVector)
    # QuadrupedLinear12 configuration
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
	error("Not Implemented")

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
