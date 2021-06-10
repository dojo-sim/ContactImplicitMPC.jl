using Random
using LinearAlgebra


################################################################################
# Test Utils
################################################################################

function get_initialization(ref_traj, t)
	Random.seed!(10)
	z = deepcopy(ref_traj.z[t])
	z += rand(length(z))
	θ = deepcopy(ref_traj.θ[t])
	return z, θ
end


################################################################################
# Test data
################################################################################
s = get_simulation("quadruped", "flat_2D_lc", "flat")
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

model = s.model
env = s.env
t = 10
nz = num_var(model, env)
nθ = num_data(model)


# # Dimensions and indices
# ix, iy1, iy2 = linearization_var_index(model, env)
# idyn, irst, ibil, ialt = linearization_term_index(model, env)
# nz = num_var(model, env)
# nθ = num_data(model)
# nx = length(ix)
# ny = length(iy1)
#
# r = zeros(nz)
# rz = zeros(nz, nz)
# rθ = zeros(nz, nθ)
#
# z = ref_traj.z[10]
# θ = ref_traj.θ[10]
# κ = 1e-6
# s.res.r!(r, z, θ, κ)
# s.res.rz!(rz, z, θ)
# s.res.rθ!(rθ, z, θ)
#
#
# r[ibil] .- (z[iy1] .* z[iy2] .- κ)
#
# # Problem data
# RZ = rz[[idyn; irst; ibil], [ix; iy1; iy2]]
# Rθ = rθ[[idyn; irst; ibil], :]
#
# E = rz[idyn, ix]
# F = rz[idyn, iy1]
# G = rz[irst, ix]
# H = rz[irst, iy1]
# J = rz[irst, iy2]
# b = r[idyn]
# c = r[irst]
#
# # Dimensions
# n = nx
# m = ny

################################################################################
# Test Interior Point on the full non linear problem
################################################################################

z1, θ1 = get_initialization(ref_traj, t)
ip1 = interior_point(z1, θ1,
    idx_ineq = inequality_indices(model, env),
    idx_soc = soc_indices(model, env),
    r! = s.res.r!,
    rm! = s.res.rm!,
    rz! = s.res.rz!,
    rθ! = s.res.rθ!,
    rz = s.rz,
    rθ = s.rθ,
    opts = InteriorPointOptions(
		r_tol=1e-8,
		κ_init=1e-8,
		κ_tol=2e-8,
		verbose=false))
@time interior_point!(ip1)
r1 = zeros(nz)
s.res.r!(r1, ip1.z, ip1.θ, 0.0)
@test norm(r1, Inf) < 1e-7

@benchmark begin
	z1, θ1 = get_initialization(ref_traj, t)
	interior_point!(ip1, z1, θ1)
end


################################################################################
# Test Mehrotra on the full non linear problem
################################################################################

z2, θ2 = get_initialization(ref_traj, t)
ip2 = mehrotra(z2, θ2,
	iy1 = linearization_var_index(model, env)[2],
	iy2 = linearization_var_index(model, env)[3],
	# ibil = linearization_term_index(model, env)[3],
    idx_ineq = inequality_indices(model, env),
    idx_soc = soc_indices(model, env),
	r! = s.res.r!,
    rm! = s.res.rm!,
    rz! = s.res.rz!,
    rθ! = s.res.rθ!,
    rz = s.rz,
    rθ = s.rθ,
    opts = Mehrotra19Options(
        max_iter_inner=100,
        r_tol=1e-8,
        κ_tol=2e-8,
		# verbose=true
		))
@time mehrotra!(ip2)
r2 = zeros(nz)
s.res.r!(r2, ip2.z, ip2.θ, 0.0)
@test norm(r2, Inf) < 1e-8

@benchmark begin
	z2, θ2 = get_initialization(ref_traj, t)
	mehrotra!(ip2, z2, θ2)
end

@profiler for k = 1:100
	z2, θ2 = get_initialization(ref_traj, t)
	mehrotra!(ip2, z2, θ2)
end



#
# ix, iy1, iy2 = linearization_var_index(model, env)
#
# @show norm(z_y1 - z[iy1])
# @show norm(z_y2 - z[iy2])
# z .+= rand(length(z))
# @show norm(z_y1 - z[iy1])
# @show norm(z_y2 - z[iy2])
#
# @show norm(Δaff_y1 - Δaff[iy1])
# @show norm(Δaff_y2 - Δaff[iy2])
# Δaff .+= rand(length(z))
# @show norm(Δaff_y1 - Δaff[iy1])
# @show norm(Δaff_y2 - Δaff[iy2])
#
# @show norm(Δ_y1 - Δ[iy1])
# @show norm(Δ_y2 - Δ[iy2])
# Δ .+= rand(length(z))
# @show norm(Δ_y1 - Δ[iy1])
# @show norm(Δ_y2 - Δ[iy2])
#
#
# w1 = z[ix]
# w2 = z[iy1]
# w3 = z[iy2]
# @show norm(w2 - z_y1)
# @show norm(w3 - z_y2)
#
# Δw1aff = Δaff[ix]
# Δw2aff = Δaff[iy1]
# Δw3aff = Δaff[iy2]
# @show norm(Δw2aff - Δaff_y1)
# @show norm(Δw3aff - Δaff_y2)
# Δw1 = Δ[ix]
# Δw2 = Δ[iy1]
# Δw3 = Δ[iy2]
# @show norm(Δw2 - Δ_y1)
# @show norm(Δw3 - Δ_y2)


################################################################################
# Test Interior Point on the linearized problem
################################################################################

im_traj1 = ImplicitTraj(ref_traj, s;
	κ = 1e-8,
	opts = InteriorPointOptions(
			κ_init = 1e-8,
			κ_tol = 2.0 * 1e-8,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver=:empty_solver,
			verbose=true))

z1, θ1 = get_initialization(ref_traj, t)
ip1 = deepcopy(im_traj1.ip[10])
@time interior_point!(ip1, z1, θ1)
r1 = zeros(nz)
s.res.r!(r1, ip1.z, ip1.θ, 0.0)
@test norm(r1, Inf) < 1e-7

@benchmark begin
	z1, θ1 = get_initialization(ref_traj, t)
	ip1 = deepcopy(im_traj1.ip[10])
	interior_point!(ip1, z1, θ1)
end

@profiler for k = 1:30000
	z1, θ1 = get_initialization(ref_traj, t)
	ip1 = deepcopy(im_traj1.ip[10])
	interior_point!(ip1, z1, θ1)
end



################################################################################
# Test Mehrotra on the linearized problem
################################################################################

im_traj2 = MehrotraImplicitTraj(ref_traj, s;
	κ = 1e-8,
	opts = Mehrotra19Options(
			κ_tol = 2.0 * 1e-8,
			r_tol = 1.0e-8,
			diff_sol = true,
			max_iter_inner=100,
			ϵ_min=0.05,
			solver=:empty_solver,
			# verbose=true
			))
z2, θ2 = get_initialization(ref_traj, t)
ip2 = deepcopy(im_traj2[10])
@time mehrotra!(ip2, z2, θ2)
r2 = zeros(nz)
s.res.r!(r2, ip2.z, ip2.θ, 0.0)
@test norm(r2, Inf) < 1e-8

@benchmark begin
	z2, θ2 = get_initialization(ref_traj, t)
	ip2 = deepcopy(im_traj2[10])
	mehrotra!(ip2, z2, θ2)
end

@profiler for k = 1:10000
	z2, θ2 = get_initialization(ref_traj, t)
	ip2 = deepcopy(im_traj2[10])
	mehrotra!(ip2, z2, θ2)
end







function MehrotraImplicitTraj(ref_traj::ContactTraj, s::Simulation;
	κ = ref_traj.κ[1],
	max_time = 60.0,
	opts = Mehrotra19Options(
			κ_init = κ[1],
			κ_tol = 2.0 * κ[1],
			r_tol = 1.0e-8,
			diff_sol = true,
			solver=:empty_solver,
			max_time=max_time,
			verbose=false))

	model = s.model
	env = s.env

	H = ref_traj.H

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nd = nq + nc + nb
	nz = num_var(model, env)
	nθ = num_data(model)

	lin = [LinearizedStep(s, ref_traj.z[t], ref_traj.θ[t], κ) for t = 1:H]

	ip =  [mehrotra(zeros(num_var(model, env)), zeros(num_data(model)),
			 idx_ineq = inequality_indices(model, env),
			 iy1 = linearization_var_index(model, env)[2],
		 	 iy2 = linearization_var_index(model, env)[3],
			 ibil = linearization_term_index(model, env)[3],
			 r! = r!,
			 rm! = rm!,
			 rz! = rz!,
			 rθ! = rθ!,
			 r   = RLin(s, lin[t].z, lin[t].θ, lin[t].r, lin[t].rz, lin[t].rθ),
			 rm  = RLin(s, lin[t].z, lin[t].θ, lin[t].r, lin[t].rz, lin[t].rθ),
			 rz = RZLin(s, lin[t].rz),
			 rθ = RθLin(s, lin[t].rθ),
			 v_pr = view(zeros(1,1), 1,1),
			 v_du = view(zeros(1,1), 1,1),
			 opts = opts) for t = 1:H]
	return ip
end
