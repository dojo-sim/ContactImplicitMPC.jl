function cable_transform(y, z)
    v1 = [0.0, 0.0, 1.0]
    # if norm(z) > norm(y)
    #     v2 = y[1:3,1] - z[1:3,1]
    # else
    #     v2 = z[1:3,1] - y[1:3,1]
    # end
    v2 = y[1:3,1] - z[1:3,1]
    normalize!(v2)
    ax = cross(v1, v2)
    ang = acos(v1'*v2)
    R = AngleAxis(ang, ax...)

    if any(isnan.(R))
        R = I
    else
        nothing
    end

    compose(Translation(z), LinearMap(R))
end

function default_background!(vis)
    setvisible!(vis["/Background"], true)
    setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
    setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
    setvisible!(vis["/Axes"], false)
end

function plot_surface!(vis::Visualizer, env::Environment{R3};  col=zeros(3), α=0.4, n::Int=50)
    f(x) = x[3] - env.surf(x[1:2])
    xlims = [-1.0, 5.0]
    ylims = [-2.0, 2.0]
    plot_surface!(vis, f, xlims=xlims, ylims=ylims, col=col, α=α, n=n)
    return nothing
end

function plot_surface!(vis::Visualizer, env::Environment{R2};  col=zeros(3), α=0.4, n::Int=50)
    f(x) = x[3] - env.surf(x[1:1])
    xlims = [-1.0, 5.0]
    ylims = [-0.1, 0.1]
    plot_surface!(vis, f, xlims=xlims, ylims=ylims, col=col, α=α, n=n)
    return nothing
end

function plot_surface!(vis::Visualizer, f::Any; xlims = [-1.0, 5.0],
        ylims = [-0.1, 0.1], col=zeros(3), α=0.4, n::Int=50)
    mesh = GeometryBasics.Mesh(f,
    	HyperRectangle(Vec(xlims[1], ylims[1], -2.0), Vec(xlims[2]-xlims[1], ylims[2]-ylims[1], 4.0)),
        MarchingCubes(), samples=(n, n, n))
    setobject!(vis["surface"], mesh,
    		   MeshPhongMaterial(color=RGBA{Float32}(col..., α)))
    return nothing
end

function get_line_material(size::Real)
    orange_col = [1,153/255,51/255]
    blue_col = [51/255,1,1]
    black_col = [0,0,0]
    orange_mat = LineBasicMaterial(color=color=RGBA(orange_col...,1.0), linewidth=size)
    blue_mat = LineBasicMaterial(color=color=RGBA(blue_col...,1.0), linewidth=size)
    black_mat = LineBasicMaterial(color=color=RGBA(black_col...,1.0), linewidth=size)
    return orange_mat, blue_mat, black_mat
end
