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

function cast3d(v::AbstractVector{T}) where {T}
	l = length(v)
	if l == 3
		return v
	elseif l == 2
		return [v[1], 0.0, v[2]]
	else
		error("Wrong vector length")
	end
end

function bivector_rotation(a::AbstractVector{T}, b::AbstractVector{T}) where {T}
	na = norm(a)
	nb = norm(b)
	v = cross(a,b)/(na*nb)
	s = norm(v)
	c = (a'*b)/(na*nb)
	vx = skew(v)
	α = (1 - c)/s^2
	R = I + vx + vx*vx*α
	return R
end

function default_background!(vis; grid::Bool=false, axes::Bool=false)
    setvisible!(vis["/Background"], true)
    setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
    setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], axes)
    setvisible!(vis["/Grid"], grid)
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

function get_material()
    orange_col = [1,153/255,51/255]
    blue_col = [51/255,1,1]
    black_col = [0,0,0]
    orange_mat = MeshPhongMaterial(color=RGBA(orange_col...,1.0))
    blue_mat = MeshPhongMaterial(color=RGBA(blue_col...,1.0))
    black_mat = MeshPhongMaterial(color=RGBA(black_col...,1.0))
    return orange_mat, blue_mat, black_mat
end
