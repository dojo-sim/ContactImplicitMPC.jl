"""
    second-order cone
"""
function second_order_cone_projection(z)
    n = length(z)

    z0 = z[1]
    z1 = view(z, 2:n)

    if norm(z1) <= z0
        return z, true
    elseif norm(z1) <= -z0
        return zero(z), false
    else
        a = 0.5 * (1.0 + z0 / norm(z1))
        z_proj = zero(z)
        z_proj[1:end-1] = a * z1
        z_proj[end] = a * norm(z1)
        return z_proj, false
    end
end

function second_order_cone_product(z, s)
    n = length(z)
    SVector{n}([z' * s; z[1] * view(s, 2:n) + s[1] * view(z, 2:n)])
end

# check that second-order cone constraints are satisfied
function second_order_cone_check(x, idx_soc)
    for idx in idx_soc
        !second_order_cone_projection(view(x, idx))[2] && return false
    end
    return true
end

# check that inequality (non-negative orthant) constraints are satisfied
function inequality_check(x, idx_ineq)
    for i in idx_ineq
        if x[i] <= 0.0
            return false
        end
    end
    return true
end

function cone_check(x, idx_ineq, idx_soc)
    !inequality_check(x, idx_ineq) && return true
    !second_order_cone_check(x, idx_soc) && return true
end


function cone_plot(λ::AbstractVector{T}, Δ::AbstractVector{T};
		ϵ = 1e-20, N::Int = 100, verbose::Bool = false, show_α::Bool = false) where {T}
	plt = plot(aspect_ratio = 1.0, fmt = :svg)

	# unit circle
	θ = Vector(0:2π/N:(1.0+1/N)*2π)
	circ = [cos.(θ) , sin.(θ)]
	plot!(plt, circ[1], circ[2], label = false)

	# Keypoints
	λp = project(λ)
	λΔp = project(λ + Δ)
	verbose && @show soc_value(λ )
	verbose && @show soc_value(λ + Δ)
	scatter!(plt, λΔp[1:1],  λΔp[2:2],  markersize = 4.0, color = :orange, label = "λ+Δ")
	scatter!(plt, λp[1:1],   λp[2:2],   markersize = 8.0, color = :red,    label = "λ")

	if show_α
		α = soc_step_length(λ, Δ, τ = 1.0, ϵ = ϵ)
		λαΔp = project(λ + α * Δ)
		verbose && @show α
		verbose && @show soc_value(λ + α * Δ)
		scatter!(plt, λαΔp[1:1], λαΔp[2:2], markersize = 2.0, color = :black,   label = "λ+αΔ")
	end

	# Line
	A = Vector(0:1/N:1.0)
	L = [[project(λ + a * Δ)[1] for a in A], [project(λ + a * Δ)[2] for a in A]]
	plot!(plt, L[1], L[2], label = false)

	display(plt)
	return nothing
end

function project(x::AbstractVector{T}; ϵ::T = 1e-15) where {T}
	if length(x) == 3
		project_3D(x, ϵ = ϵ)
	elseif length(x) == 2
		project_2D(x, ϵ = ϵ)
	end
end

function project_3D(x::AbstractVector{T}; ϵ::T = 1e-15) where {T}
	@assert length(x) == 3
	xs = x[1]
	xv = x[2:end]
	r = norm(xv) / (xs + ϵ)
	v = xv ./ norm(xv)
	p = r * v
	return p
end

function project_2D(x::AbstractVector{T}; ϵ::T = 1e-15) where {T}
	@assert length(x) == 2
	xs = x[1]
	xv = x[2:end]
	r = norm(xv) / (xs + ϵ)
	v = xv ./ norm(xv)
	p = r * v
	return [p[1], 0.0]
end
