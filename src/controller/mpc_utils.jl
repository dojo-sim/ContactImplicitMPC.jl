function rot_n_stride!(traj::ContactTraj, stride::SizedArray)
    rotate!(traj)
    mpc_stride!(traj, stride)
    return nothing
end

function rotate!(traj::ContactTraj)
    H = traj.H
	# TODO: remove copy
    q1 = copy(traj.q[1])
    u1 = copy(traj.u[1])
    w1 = copy(traj.w[1])
    γ1 = copy(traj.γ[1])
    b1 = copy(traj.b[1])
    z1 = copy(traj.z[1])
    θ1 = copy(traj.θ[1])

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

    traj.q[end] .= q1
    traj.u[end] .= u1
    traj.w[end] .= w1
    traj.γ[end] .= γ1
    traj.b[end] .= b1
    traj.z[end] .= z1
    traj.θ[end] .= θ1

    return nothing
end

"""
    Update the last two cofiguaraton to be equal to the fisrt two up to an constant offset.
"""
function mpc_stride!(traj::ContactTraj, stride::SizedArray)
    H = traj.H

    for t = H+1:H+2
        traj.q[t] .= deepcopy(traj.q[t-H] + stride)
        update_z!(traj, t - 2)
        update_θ!(traj, t - 2)
        update_θ!(traj, t - 3)
    end

    return nothing
end

function get_stride(model::ContactDynamicsModel, traj::ContactTraj) #TODO: dispatch over environment / model
    stride = zeros(SizedVector{model.dim.q})
    stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return stride
end

# TODO: make more efficient / allocation free
function update_altitude!(alt, model::ContactDynamicsModel, traj, t, N_sample;
	threshold = 1.0, verbose = false)

	# idx1 = max(0, t - N_sample) + 1
	idx1 = max(0, t-1 - N_sample) + 1

	for i = 1:model.dim.c
		γ_max = 0.0
		idx_max = 0

		# for j = idx1:t
		for j = idx1:t-1
			if traj.γ[j][i] > γ_max
				γ_max =  traj.γ[j][i]
				idx_max = j
			end
		end
		# @show ([γ[1] for γ in traj.γ[max(1, t-3): t]])
		if γ_max > threshold
			alt[i] = ϕ_func(model, traj.q[idx_max+2])[i] #TODO check this this is correct
			# alt[i] = ϕ_func(model, traj.q[idx_max+2])[i]
			verbose && println(" ")
			verbose && println("point $i in contact")
			verbose && println("sim_step : $idx_max")
			verbose && println("alt      : $(ϕ_func(model, traj.q[idx_max])[i])")
			verbose && println("force    : $(traj.γ[idx_max][i])")
		end
	end
end

function live_plotting(model::ContactDynamicsModel, ref_traj::ContactTraj,
		sim_traj::ContactTraj, newton::Newton, q0::AbstractVector{T},
		q1::AbstractVector{T}, t::Int) where {T}
	nq = model.dim.q
	nu = model.dim.u

	ql = 1
	qu = nq
	ul = 1
	uu = nu
	plt = plot(layout=grid(2,1,heights=[0.7, 0.3], figsize=[(1000, 1000),(400,400)]), legend=false, xlims=(0,20))
	plot!(plt[1,1], hcat([Vector(x[ql:qu]) for x in newton.traj.q]...)', color=:blue, linewidth=1.0)
	plot!(plt[1,1], hcat([Vector(x[ql:qu]) for x in ref_traj.q]...)', linestyle=:dot, color=:red, linewidth=3.0)

	scatter!((2-1/N_sample)ones(nq), sim_traj.q[t+0], markersize=8.0, color=:lightgreen)
	scatter!(2*ones(nq), sim_traj.q[t+1], markersize=8.0, color=:lightgreen)
	# scatter!(plt[1,1], 3*ones(nq), sim_traj.q[t+2], markersize=8.0, color=:lightgreen)

	# @show norm(sim_traj.q[t+1] - q1, Inf)
	scatter!(plt[1,1], 1*ones(nq), q0, markersize=6.0, color=:blue)
	scatter!(plt[1,1], 2*ones(nq), q1, markersize=6.0, color=:blue)

	scatter!(plt[1,1], 1*ones(nq), ref_traj.q[1], markersize=4.0, color=:red)
	scatter!(plt[1,1], 2*ones(nq), ref_traj.q[2], markersize=4.0, color=:red)
	scatter!(plt[1,1], 3*ones(nq), ref_traj.q[3], markersize=4.0, color=:red)

	plot!(plt[2,1], hcat([Vector(x[ul:uu]) for x in newton.traj.u]...)', color=:blue, linewidth=1.0)
	plot!(plt[2,1], hcat([Vector(x[ul:uu]) for x in ref_traj.u[1:end]]...)', linestyle=:dot, color=:red, linewidth=3.0)
	display(plt)
end
