include(joinpath(pwd(),"src/dynamics/flamingo/visuals.jl"))
vis = Visualizer()
render(vis)

dir = joinpath(pwd(), "src/dynamics/flamingo")
urdf = joinpath(dir, "mesh", "flamingo.urdf")
mechanism = MeshCatMechanisms.parse_urdf(urdf, remove_fixed_tree_joints = true)
mechanism_visuals = URDFVisuals(urdf, package_path=[dir])
mvis = MechanismVisualizer(mechanism, mechanism_visuals, vis[:Flamingo])

function set_alpha!(visuals::Vector{VisualElement}, α)
    for el in visuals
        c = el.color
        c_new = RGBA(red(c),green(c),blue(c),α)
        el.color = c_new
    end
end

vis_el = visual_elements(mechanism, mechanism_visuals)
set_alpha!(vis_el, 0.25)

mvis2 = MechanismVisualizer(mechanism, vis[Symbol("shadow")])
MeshCatMechanisms._set_mechanism!(mvis2, vis_el)
MeshCatMechanisms._render_state!(mvis2)

setvisible!(vis[:Flamingo], false)
# setvisible!(vis[:ghost], false)
setvisible!(vis[:shadow], true)
