# using Colors
# using CoordinateTransformations
# using FileIO
# using GeometryBasics
# using MeshCat, MeshIO, Meshing
# using Rotations
#
function visualize!(vis, model::Particle, q;
	Δt = 0.1, r = 0.25,
	name = "1",
	color = RGBA(1.0, 0.0, 0.0, 1.0))

	default_background!(vis)

    setobject!(vis["particle" * name],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = color))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["particle" * name], MeshCat.Translation(q[t][1:3]...))
        end
    end

	settransform!(vis["/Cameras/default"],
	    compose(Translation(-2.5, 7.5, 1.0),LinearMap(RotZ(0.0))))

    MeshCat.setanimation!(vis, anim)
end
#
#




using Colors
using CoordinateTransformations
using FileIO
using GeometryBasics
using MeshCat, MeshIO, Meshing
using Rotations
#
# function visualize!(vis, model::Particle, q;
# 	Δt = 0.1, r = 0.25)
#
# 	default_background!(vis)
#
#     setobject!(vis["particle"],
# 		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
# 		convert(Float32, r)),
# 		MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, 1.0)))
#
#     anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))
#
#     for t = 1:length(q)
#         MeshCat.atframe(anim, t) do
#             settransform!(vis["particle"], MeshCat.Translation(q[t][1], q[t][2], q[t][3]))
#         end
#     end
#
#     MeshCat.setanimation!(vis, anim)
# end
#


function plot_lines!(vis::Visualizer, model::Particle, q::AbstractVector;
		r=0.25, size=10, name::Symbol = :particle, col::Bool=true)
	# Point Traj
	com_point = Vector{Point{3,Float64}}()

	for qi in q
		push!(com_point, Point(qi[1], qi[2], qi[3] + r))
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

function build_robot!(vis::Visualizer, model::Particle; name::Symbol = :Particle, r=0.25, α=1.0)
	nc = model.dim.c
	r = convert(Float32, r)
	body_mat = MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, α))

	default_background!(vis)

    setobject!(vis[name][:robot]["particle"], GeometryBasics.Rect(
		-r * Vec(1,1,1.),
		2r * Vec(1,1,1.)),
		body_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::Particle, q::AbstractVector; name::Symbol=:Particle, r=0.25)
	nc = model.dim.c
	settransform!(vis[name][:robot]["particle"], Translation(q[1], q[2], q[3]+r))
	return nothing
end

function contact_point(model::Particle, q::AbstractVector; r=0.25)
	return [q[1:3]]
end
