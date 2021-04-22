using MeshCat
using MeshCatMechanisms
vis = Visualizer()
# open(vis)
render(vis)
# default_background!(vis)

flamingo = MeshCatMechanisms.URDFVisuals(joinpath(@__DIR__, "mesh", "flamingo.urdf"))

urdf = joinpath(@__DIR__, "mesh", "flamingo.urdf")
mechanism = MeshCatMechanisms.parse_urdf(urdf)
mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf, package_path=[@__DIR__]), vis[:Flamingo])
#
# cube = HyperRectangle(-.1725, +0.12, -0.87, 0.1725, 0.12, +0.42)
# mat = MeshPhongMaterial(color=RGBA(1.0,0.0,0.0,1.0))
# setobject!(vis[:cube], cube, mat)
#
# model = flamingo
#
# q = [0.0, 0.83, 0.05, 0.3, -0.3, 0.1, -0.5, -π/2, -π/2]
# build_robot!(vis, model)
# set_robot!(vis, model, q)
