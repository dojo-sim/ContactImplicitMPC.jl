using Colors
using CoordinateTransformations
using FileIO
using GeometryBasics
using MeshCat, MeshIO, Meshing
using Rotations

function visualize!(vis, model::Particle2D, q;
	Δt = 0.1, r = 0.25)

	default_background!(vis)

    setobject!(vis["particle"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["particle"], MeshCat.Translation(q[t][1], 0.0, q[t][2]))
        end
    end

    MeshCat.setanimation!(vis, anim)
end



function plot_lines!(vis::Visualizer, model::Particle2D, q::AbstractVector;
		r=0.25, size=10, name::Symbol=:particle_2D, col::Bool=true)
	# Point Traj
	com_point = Vector{Point{3,Float64}}()

	for qi in q
		push!(com_point, Point(qi[1], 0.0, qi[2] + r))
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["lines"]["com"], MeshCat.Line(com_point, orange_mat))
	else
		setobject!(vis[name]["lines"]["com"], MeshCat.Line(com_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Particle2D; name::Symbol=:Particle2D, r=0.25, α=1.0)
	nc = model.dim.c
	r = convert(Float32, r)
	body_mat = MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, α))

	default_background!(vis)

    setobject!(vis[name][:robot]["particle_2D"], GeometryBasics.Rect(
		Vec(-1.0 * r, -1.0 * r, -1.0 * r),
		Vec(2.0 * r, 2.0 * r, 2.0 * r)),
		body_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::Particle2D, q::AbstractVector; name::Symbol=:Particle2D, r=0.25)
	nc = model.dim.c

	settransform!(vis[name][:robot]["particle_2D"], Translation(q[1], 0.0, q[2]+r))
	return nothing
end

function contact_point(model::Particle2D, q::AbstractVector; r=0.25)
	pc = [q[1], 0.0, q[2]]
	return [pc]
end
