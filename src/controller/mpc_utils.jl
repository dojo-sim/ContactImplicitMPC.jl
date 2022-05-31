function rot_n_stride!(traj::ContactTraj, cache::ContactTraj, stride::Vector{T}, window::Vector{Int}) where T
    rotate!(traj, cache)
    mpc_stride!(traj, stride, window)
    return nothing
end

function rotate!(traj::ContactTraj, cache::ContactTraj)
    H = traj.H

	cache.q[1] .= traj.q[1]
    cache.u[1] .= traj.u[1]
    cache.w[1] .= traj.w[1]
    cache.γ[1] .= traj.γ[1]
    cache.b[1] .= traj.b[1]
    cache.z[1] .= traj.z[1]
    cache.θ[1] .= traj.θ[1]

    for t = 1:H+1
        traj.q[t] .= traj.q[t+1]
    end

    for t = 1:H-1
        traj.u[t] .= traj.u[t+1]
        traj.w[t] .= traj.w[t+1]
        traj.γ[t] .= traj.γ[t+1]
        traj.b[t] .= traj.b[t+1]
        traj.z[t] .= traj.z[t+1]
        traj.θ[t] .= traj.θ[t+1]
    end

    traj.q[end] .= cache.q[1]
    traj.u[end] .= cache.u[1]
    traj.w[end] .= cache.w[1]
    traj.γ[end] .= cache.γ[1]
    traj.b[end] .= cache.b[1]
    traj.z[end] .= cache.z[1]
    traj.θ[end] .= cache.θ[1]

    return nothing
end

function set_window!(traj::ContactTraj, ref::ContactTraj, window)
    
	for (i, t) in enumerate(window)
		traj.q[i] .= ref.q[t]
	end

	for (i, t) in enumerate(window[1:end-2])
		traj.u[i] .= ref.u[t]
		traj.w[i] .= ref.w[t]
		traj.γ[i] .= ref.γ[t]
		traj.b[i] .= ref.b[t]
		traj.z[i] .= ref.z[t]
		traj.θ[i] .= ref.θ[t]
	end

    return nothing
end

function set_trajectory!(traj::ContactTraj, ref_traj::ContactTraj)
	H = ref_traj.H
	for t = 1:H
		traj.q[t] .= ref_traj.q[t]
		traj.u[t] .= ref_traj.u[t]
		traj.w[t] .= ref_traj.w[t]
		traj.γ[t] .= ref_traj.γ[t]
		traj.b[t] .= ref_traj.b[t]
		traj.z[t] .= ref_traj.z[t]
		traj.θ[t] .= ref_traj.θ[t]
	end
	traj.q[H+1] .= ref_traj.q[H+1]
	traj.q[H+2] .= ref_traj.q[H+2]
	return nothing
end

"""
    Update the last two cofiguaraton to be equal to the fisrt two up to an constant offset.
"""
function mpc_stride!(traj::ContactTraj, stride::Vector{T}, window::Vector{Int}) where T
    H = traj.H

    for t = H+1:H+2
		traj.q[t] .= traj.q[t-H]
		traj.q[t] .+= stride

		# update_z!(traj, t - 2)
		τ = t - 2
		traj.z[τ][traj.iq2] = traj.q[τ+2]

		# update_θ!(traj, t - 2)
		traj.θ[τ][traj.iq0] = traj.q[τ]
		traj.θ[τ][traj.iq1] = traj.q[τ+1]

		# update_θ!(traj, t - 3)
		τ = t - 3
		traj.θ[τ][traj.iq0] = traj.q[τ]
		traj.θ[τ][traj.iq1] = traj.q[τ+1]
    end

    return nothing
end

function get_stride(model::Model, traj::ContactTraj) #TODO: dispatch over environment / model
    stride = zeros(model.nq)
    stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return stride
end

function update_altitude!(alt::Vector{T}, ϕ::Vector{T}, s::Simulation{T}, traj::Trajectory{T}, t::Int, nc::Int, N_sample::Int;
	threshold = 1.0, verbose = false) where T

	idx1 = max(0, t - 1 - N_sample) + 1

	for i = 1:nc
		γ_max = 0.0
		idx_max = 0

		for j = idx1:t-1
			if traj.γ[j][i] > γ_max
				γ_max = traj.γ[j][i]
				idx_max = j
			end
		end

		if γ_max > threshold
			s.ϕ(ϕ, traj.q[idx_max+2])
			alt[i] = ϕ[i]
			verbose && println(" ")
			verbose && println("point $i in contact")
			verbose && println("sim_step : $idx_max")
			verbose && println("alt      : $ϕ")
			verbose && println("force    : $(traj.γ[idx_max][i])")
		end
	end
end

function update_altitude!(alt::Vector{T}, ϕ::Vector{T}, s::Simulation{T}, x::Vector{T}, nq::Int, nc::Int;
	threshold = 1.0, 
	verbose = false) where T

	s.ϕ(ϕ, x[1:2nq])
	γ = x[2nq .+ (1:nc)]

	for i = 1:nc
		if γ[i] > threshold
			alt[i] = ϕ[i]
			verbose && println(" ")
			verbose && println("point $i in contact")
			verbose && println("sim_step : $idx_max")
			verbose && println("alt      : $(ϕ[i])")
			verbose && println("force    : $(γ[i])")
		end
	end
end

function live_plotting(model::Model{T}, ref_traj::ContactTraj,
		sim_traj::Trajectory, newton::Newton, q0::AbstractVector{T},
		q1::AbstractVector{T}, t::Int) where T
	nq = model.nq
	nu = model.nu

	ql = 1
	qu = nq
	ul = 1
	uu = nu
	plt = plot(layout=grid(2,1,heights=[0.7, 0.3], figsize=[(1000, 1000),(400,400)]), legend=:bottomright, xlims=(0,20))
	plot!(plt[1,1], hcat([Vector(x[ql:qu]) for x in newton.traj.q]...)', color=:blue, linewidth=1.0, label=nothing)
	plot!(plt[1,1], hcat([Vector(x[ql:qu]) for x in ref_traj.q]...)', linestyle=:dot, color=:red, linewidth=3.0, label=nothing)

	scatter!((2-1/N_sample)ones(nq), sim_traj.q[t+0], markersize=8.0, color=:lightgreen, label="simul. q0 q1")
	scatter!(2*ones(nq), sim_traj.q[t+1], markersize=8.0, color=:lightgreen, label=nothing)

	scatter!(plt[1,1], 1*ones(nq), q0, markersize=6.0, color=:blue, label="newton q0 q1")
	scatter!(plt[1,1], 2*ones(nq), q1, markersize=6.0, color=:blue, label=nothing)

	scatter!(plt[1,1], 1*ones(nq), ref_traj.q[1], markersize=4.0, color=:red, label="refer. q0 q1 q2")
	scatter!(plt[1,1], 2*ones(nq), ref_traj.q[2], markersize=4.0, color=:red, label=nothing)
	scatter!(plt[1,1], 3*ones(nq), ref_traj.q[3], markersize=4.0, color=:red, label=nothing)

	plot!(plt[2,1], hcat([Vector(x[ul:uu]) for x in newton.traj.u]...)', color=:blue, linewidth=1.0, label=nothing)#"u newton")
	plot!(plt[2,1], hcat([Vector(x[ul:uu]) for x in ref_traj.u[1:end]]...)', linestyle=:dot, color=:red, linewidth=3.0, label=nothing)#"u refer.")
	display(plt)
end

