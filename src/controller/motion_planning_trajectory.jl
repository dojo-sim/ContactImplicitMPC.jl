struct Trajectory
	x

	τ
	ν

	q
	u
	w
	γ
	b

	z
	θ

	q0
	q1

	zq2
	zγ1
	zb1
	zψ1
	zη1
	zs1
	zs2

	θq0
	θq1
	θu1
	θw1

	H
	h
	κ
end

function trajectory(model::ContactDynamicsModel, q0, q1, H, h, κ,
	w = [@SVector zeros(model.dim.w) for t = 1:H])

	dim = model.dim
	nq = dim.q
    nu = dim.u
    nw = dim.w
    nc = dim.c
    nb = dim.b
	nz = num_var(dim)
	nθ = num_data(dim)
	nr = nq + nu + nc + nb
	nd = nq + nc + nb
	nτ = nr * H
	nν = nd * H

	x = zeros(nτ + nν)
	τ = [view(x, (t - 1) * nr .+ (1:nr)) for t = 1:H]
	ν = [view(x, nτ + (t - 1) * nd .+ (1:nd)) for t = 1:H]

	q = [t == 1 ? view(q0, 1:nq) : (t == 2 ? view(q1, 1:nq) : view(τ[t-2], 1:nq)) for t = 1:H+2]
	u = [view(τ[t], nq .+ (1:nu)) for t = 1:H]
	w = [zeros(nw) for k = 1:H]
	γ = [view(τ[t], nq + nu .+ (1:nc)) for t = 1:H]
	b = [view(τ[t], nq + nu + nc .+ (1:nb)) for t = 1:H]

	z = [zeros(nz) for t = 1:H]
	θ = [zeros(nθ) for t = 1:H]

	zq2 = [view(z[t], 1:nq) for t = 1:H]
	zγ1 = [view(z[t], nq .+ (1:nc)) for t = 1:H]
	zb1 = [view(z[t], nq + nc .+ (1:nb)) for t = 1:H]
	zψ1 = [view(z[t], nq + nc + nb .+ (1:nc)) for t = 1:H]
	zη1 = [view(z[t], nq + nc + nb + nc .+ (1:nb)) for t = 1:H]
	zs1 = [view(z[t], nq + nc + nb + nc + nb .+ (1:nc)) for t = 1:H]
	zs2 = [view(z[t], nq + nc + nb + nc + nb + nc .+ (1:nc)) for t = 1:H]

	θq0 = [view(θ[t], 1:nq) for t = 1:H]
	θq1 = [view(θ[t], nq .+ (1:nq)) for t = 1:H]
	θu1 = [view(θ[t], 2nq .+ (1:nu)) for t = 1:H]
	θw1 = [view(θ[t], 2nq + nu .+ (1:nw)) for t = 1:H]

	Trajectory(x, τ, ν,	q, u, w, γ, b, z, θ, q0, q1, zq2, zγ1, zb1, zψ1, zη1, zs1, zs2,
		θq0, θq1, θu1, θw1,	H, h, κ)
end

function trajectory_x(model::ContactDynamicsModel, x, q0, q1, H, h, κ,
	w = [@SVector zeros(model.dim.w) for t = 1:H])

	dim = model.dim
	nq = dim.q
    nu = dim.u
    nw = dim.w
    nc = dim.c
    nb = dim.b
	nz = num_var(dim)
	nθ = num_data(dim)
	nr = nq + nu + nc + nb
	nd = nq + nc + nb
	nτ = nr * H
	nν = nd * H

	τ = [view(x, (t - 1) * nr .+ (1:nr)) for t = 1:H]
	ν = [view(x, nτ + (t - 1) * nd .+ (1:nd)) for t = 1:H]

	q = [t == 1 ? view(q0, 1:nq) : (t == 2 ? view(q1, 1:nq) : view(τ[t-2], 1:nq)) for t = 1:H+2]
	u = [view(τ[t], nq .+ (1:nu)) for t = 1:H]
	w = [zeros(nw) for k = 1:H]
	γ = [view(τ[t], nq + nu .+ (1:nc)) for t = 1:H]
	b = [view(τ[t], nq + nu + nc .+ (1:nb)) for t = 1:H]

	z = [zeros(nz) for t = 1:H]
	θ = [zeros(nθ) for t = 1:H]

	zq2 = [view(z[t], 1:nq) for t = 1:H]
	zγ1 = [view(z[t], nq .+ (1:nc)) for t = 1:H]
	zb1 = [view(z[t], nq + nc .+ (1:nb)) for t = 1:H]
	zψ1 = [view(z[t], nq + nc + nb .+ (1:nc)) for t = 1:H]
	zη1 = [view(z[t], nq + nc + nb + nc .+ (1:nb)) for t = 1:H]
	zs1 = [view(z[t], nq + nc + nb + nc + nb .+ (1:nc)) for t = 1:H]
	zs2 = [view(z[t], nq + nc + nb + nc + nb + nc .+ (1:nc)) for t = 1:H]

	θq0 = [view(θ[t], 1:nq) for t = 1:H]
	θq1 = [view(θ[t], nq .+ (1:nq)) for t = 1:H]
	θu1 = [view(θ[t], 2nq .+ (1:nu)) for t = 1:H]
	θw1 = [view(θ[t], 2nq + nu .+ (1:nw)) for t = 1:H]

	Trajectory(x, τ, ν,	q, u, w, γ, b, z, θ, q0, q1, zq2, zγ1, zb1, zψ1, zη1, zs1, zs2,
		θq0, θq1, θu1, θw1,	H, h, κ)
end

function update_θ!(traj::Trajectory)
	for t = 1:traj.H
		traj.θq0[t] .= traj.q[t]
		traj.θq1[t] .= traj.q[t+1]
		traj.θu1[t] .= traj.u[t]
		traj.θw1[t] .= traj.w[t]
		traj.θ[t][end] = traj.h

		# traj.θ[t] .= [traj.q[t]; traj.q[t+1]; traj.u[t]; traj.w[t]; traj.h]
	end
	return nothing
end

function update_z!(traj::Trajectory)
	for t = 1:traj.H
		traj.zq2[t] .= traj.q[t+2]
		traj.zγ1[t] .= traj.γ[t]
		traj.zb1[t] .= traj.b[t]
		# traj.zψ1[t]
		# traj.zη1[t]
		# traj.zs1[t]
		# traj.zs2[t]
	end
	return nothing
end
