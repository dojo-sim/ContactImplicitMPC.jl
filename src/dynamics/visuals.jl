function plot_surface!(vis::Visualizer, env::Environment{R3};  col=zeros(3), α=0.4,
		n::Int=50, xlims = [-1.0, 5.0], ylims = [-2.0, 2.0])
    f(x) = x[3] - env.surf(x[1:2])[1]
    plot_surface!(vis, f, xlims=xlims, ylims=ylims, col=col, α=α, n=n)
    return nothing
end

function plot_surface!(vis::Visualizer, env::Environment{R2};  col=zeros(3), α=0.4,
		n::Int=50, xlims = [-1.0, 5.0], ylims = [-0.1, 0.1])
    f(x) = x[3] - env.surf(x[1:1])[1]
    plot_surface!(vis, f, xlims=xlims, ylims=ylims, col=col, α=α, n=n)
    return nothing
end

function plot_surface!(vis::Visualizer, f::Any; xlims = [-1.0, 5.0],
        ylims = [-0.1, 0.1], col=zeros(3), α=0.4, n::Int=200)
    mesh = GeometryBasics.Mesh(f,
    	HyperRectangle(Vec(xlims[1], ylims[1], -2.0), Vec(xlims[2]-xlims[1], ylims[2]-ylims[1], 4.0)),
        Meshing.MarchingCubes(), samples=(n, n, Int(floor(n/8))))
    setobject!(vis["surface"], mesh,
    		   MeshPhongMaterial(color=RGBA{Float32}(col..., α)))
    return nothing
end

function animate_robot!(vis::Visualizer, anim::MeshCat.Animation, model::ContactModel,
		q::AbstractVector; name::Symbol=model_name(model))
	for t in 1:length(q)
		MeshCat.atframe(anim, t) do
			set_robot!(vis, model, q[t], name=name)
		end
	end
	setanimation!(vis, anim)
	return nothing
end

function visualize_robot!(vis::Visualizer, model::ContactModel, q::AbstractVector;
		h=0.01, α=1.0,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	build_robot!(vis, model, name=name, α=α)
	animate_robot!(vis, anim, model, q, name=name)
	return anim
end

function visualize_robot!(vis::Visualizer, model::ContactModel, traj::ContactTraj;
		sample=max(1, Int(floor(traj.H / 100))), h=traj.h*sample,  α=1.0,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	visualize_robot!(vis, model, traj.q[3:sample:end]; anim=anim, name=name, h=h, α=α)
	return anim
end

function build_meshrobot!(vis::Visualizer, model::ContactModel, urdf::String,
		package_path::String; name::Symbol=model_name(model), α=1.0)
	default_background!(vis)

	mechanism = MeshCatMechanisms.parse_urdf(urdf)
    visuals = URDFVisuals(urdf, package_path=[package_path])
    state = MeshCatMechanisms.MechanismState(mechanism)
    vis_el = MeshCatMechanisms.visual_elements(mechanism, visuals)
    set_alpha!(vis_el,α)

    mvis = MechanismVisualizer(state, vis[name, :world])
    MeshCatMechanisms._set_mechanism!(mvis, vis_el)
    MeshCatMechanisms._render_state!(mvis)
	return mvis
end

function set_alpha!(visuals::Vector{VisualElement}, α)
    for el in visuals
        c = el.color
        c_new = RGBA(red(c),green(c),blue(c),α)
        el.color = c_new
    end
end

function set_meshrobot!(vis::Visualizer, mvis::MechanismVisualizer, model::ContactModel,
		q::AbstractVector; name::Symbol=model_name(model))

	q_mesh, tform = convert_config(model, q)
	set_configuration!(mvis, q_mesh)
	settransform!(vis[name, :world], tform)

	return nothing
end

function animate_meshrobot!(vis::Visualizer, mvis::MechanismVisualizer, anim::MeshCat.Animation,
		model::ContactModel, q::AbstractVector; name::Symbol=model_name(model))
	for t in 1:length(q)
		MeshCat.atframe(anim, t) do
			set_meshrobot!(vis, mvis, model, q[t], name=name)
		end
	end
	setanimation!(vis, anim)
	return nothing
end

function visualize_meshrobot!(vis::Visualizer, model::ContactModel, q::AbstractVector;
		h=0.01,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	mvis = build_meshrobot!(vis, model, name=name)
	animate_meshrobot!(vis, mvis, anim, model, q, name=name)

	return anim
end

function visualize_meshrobot!(vis::Visualizer, model::ContactModel, traj::ContactTraj;
		sample=max(1, Int(floor(traj.H / 100))), h=traj.h*sample,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	anim = visualize_meshrobot!(vis, model, traj.q[3:sample:end]; anim=anim, name=name, h=h)
	return anim
end

function build_force!(vis::Visualizer, model::ContactModel; name::Symbol=model_name(model))
	default_background!(vis)
	orange_mat, blue_mat, black_mat = get_material()
	nc = model.dim.c

	for i = 1:nc
		fri_vis = ArrowVisualizer(vis[name][:force][:friction]["$i"])
		imp_vis = ArrowVisualizer(vis[name][:force][:impact]["$i"])
		con_vis = ArrowVisualizer(vis[name][:force][:contact]["$i"])
		setobject!(imp_vis, orange_mat)
		setobject!(fri_vis, black_mat)
		setobject!(con_vis, blue_mat)
	end
	return nothing
end

function contact_force(model::ContactModel, pc::AbstractVector,
		γ::AbstractVector, b::AbstractVector;
		E::AbstractMatrix=Matrix(friction_mapping(model.env)))
	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb/nc)

	r = [Matrix(rotation(model.env, pc[i])') for i=1:nc]
	bc = [r[i] * [E*b[(i-1) * nf .+ (1:nf)];  0.0] for i = 1:nc] # TODO: make efficient
	λc = [r[i] * [E*b[(i-1) * nf .+ (1:nf)]; γ[i]] for i = 1:nc] # TODO: make efficient
	return bc, λc
end

function set_force!(vis::Visualizer, model::ContactModel, q::AbstractVector,
		γ::AbstractVector, b::AbstractVector; name::Symbol=model_name(model), shift=-0.04,
		E::AbstractMatrix=Matrix(friction_mapping(model.env)))

	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb/nc)
	p_shift = [0, shift, 0]

	pc = contact_point(model, q)
	bc, λc = contact_force(model, pc, γ, b, E=E)

	imp_vis = [ArrowVisualizer(vis[name][:force][:impact]["$i"]) for i=1:nc]
	fri_vis = [ArrowVisualizer(vis[name][:force][:friction]["$i"]) for i=1:nc]
	con_vis = [ArrowVisualizer(vis[name][:force][:contact]["$i"]) for i=1:nc]

	for i = 1:nc
		bi = cast3d(bc[i])
		λi = cast3d(λc[i])

		settransform!(imp_vis[i],
			Point(0.0, 0.0, 0.0),
			Vec(γ[i], 0.0, 0.0),
			shaft_radius=0.008,
			max_head_radius=0.020)
		settransform!(fri_vis[i],
			Point(0.0, 0.0, 0.0),
			Vec(norm(bi), 0.0, 0.0),
			shaft_radius=0.008,
			max_head_radius=0.020)
		settransform!(con_vis[i],
			Point(0.0, 0.0, 0.0),
			Vec(norm(λi), 0.0, 0.0),
			shaft_radius=0.008,
			max_head_radius=0.020)

		trans = Translation(pc[i] + p_shift)
		rot = LinearMap(rotation_3d(model.env, pc[i])')
		transform = compose(trans, rot)
		settransform!(vis[name][:force][:impact]["$i"], transform)

		rot_fri = LinearMap(bivector_rotation([0.0, 0.0, 1.0], bi))
		transform_fri = compose(trans, rot_fri)
		settransform!(vis[name][:force][:friction]["$i"], transform_fri)

		rot_con = LinearMap(bivector_rotation([0.0, 0.0, 1.0], λi))
		transform_con = compose(trans, rot_con)
		settransform!(vis[name][:force][:contact]["$i"], transform_con)
	end
	return nothing
end

function animate_force!(vis::Visualizer, anim::MeshCat.Animation, model::ContactModel,
	q::AbstractVector, γ::AbstractVector, b::AbstractVector;
	name::Symbol=model_name(model))
	E = Matrix(friction_mapping(model.env))
	for t in 1:length(q)
		MeshCat.atframe(anim, t) do
			set_force!(vis, model, q[t], γ[t], b[t], name=name, E=E)
		end
	end
	setanimation!(vis, anim)
	return nothing
end

function visualize_force!(vis::Visualizer, model::ContactModel,
		q::AbstractVector, γ::AbstractVector, b::AbstractVector; h=0.01,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	build_force!(vis, model, name=name)
	animate_force!(vis, anim, model, q, γ, b, name=name)
	return anim
end

function visualize_force!(vis::Visualizer, model::ContactModel,
		traj::ContactTraj; sample=max(1,Int(floor(traj.H/100))), h=traj.h*sample,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	visualize_force!(vis, model, traj.q[3:sample:end], traj.γ[1:sample:end], traj.b[1:sample:end];
			anim=anim, name=name, h=h)
	return anim
end
