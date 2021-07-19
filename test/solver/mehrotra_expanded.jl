using Random
using LinearAlgebra
# const ContactControl = Main

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

function interior_point_timing(ref_traj, t, ip0)
	e1 = @belapsed begin
		z, θ = get_initialization(ref_traj, t)
		ip = deepcopy(ip1)
	end
	e2 = @belapsed begin
		z, θ = get_initialization(ref_traj, t)
		ip = deepcopy(ip1)
		interior_point_solve!(ip, z, θ)
	end
	@show e1
	@show e2
	return e2 - e1
end

################################################################################
# Test data
################################################################################

s = ContactControl.get_simulation("quadruped", "flat_2D_lc", "flat")
ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(ContactControl.module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

model = s.model
env = s.env
t = 10
nz = num_var(model, env)
nθ = num_data(model)

################################################################################
# Mehrotra12
# - linear       NO
# - euclidean    YES
# - conic        NO
################################################################################

z2, θ2 = get_initialization(ref_traj, t)
ip2 = mehrotra(z2, θ2,
	ix = linearization_var_index(model, env)[1],
	iy1 = linearization_var_index(model, env)[2],
	iy2 = linearization_var_index(model, env)[3],
	idyn = linearization_term_index(model, env)[1],
	irst = linearization_term_index(model, env)[2],
	ibil = linearization_term_index(model, env)[3],
    idx_ineq = inequality_indices(model, env),
    idx_soc = soc_indices(model, env),
	r! = s.res.r!,
    rm! = s.res.rm!,
    rz! = s.res.rz!,
    rθ! = s.res.rθ!,
    rz = s.rz,
    rθ = s.rθ,
    opts = Mehrotra12Options(
        max_iter_inner=30,
        r_tol=1e-8,
        κ_tol=1e-8,
		# verbose=true
		))
interior_point_solve!(ip2, z2, θ2)
r2 = zeros(nz)
ip2.methods.r!(r2, ip2.z, ip2.θ, 0.0)
@testset "Mehrotra12 Nonlinear" begin
	@test norm(r2, Inf) < 1e-8
	@test abs(norm(r2, Inf) - 7.84e-9) < 1e-11
	@test ip2.iterations == 9
end


################################################################################
# Mehrotra12
# - linear       YES
# - euclidean    YES
# - conic        NO
################################################################################

im_traj2 = ImplicitTraj(ref_traj, s;
	κ = 1e-8,
	ip_type = :mehrotra,
	opts = Mehrotra12Options(
			κ_tol = 2.0 * 1e-8,
			r_tol = 1.0e-8,
			diff_sol = true,
			max_iter_inner=100,
			ϵ_min=0.05,
			solver=:empty_solver,
			# verbose=true
			))
z2, θ2 = get_initialization(ref_traj, t)
ip2 = deepcopy(im_traj2.ip[10])
interior_point_solve!(ip2, z2, θ2)
r2 = zeros(nz)
s.res.r!(r2, ip2.z, ip2.θ, 0.0)
@testset "Mehrotra12 Linear" begin
	@test norm(r2, Inf) < 1e-8
	@test (norm(r2, Inf) - 4.4e-9) < 1e-13
	@test ip2.iterations == 3
end




#
# ################################################################################
# # Test data
# ################################################################################
#
# s = ContactControl.get_simulation("box", "flat_3D_lc", "flat_lc")
# model = s.model
# env = s.env
#
# h = 0.05
# H = 100
# q0 = SVector{model.dim.q}([2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0])
# q1 = SVector{model.dim.q}([2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0])
#
# sim = simulator(s, q0, q1, h, H,
# 	space = rquat_space,
#     ip_opts = ContactControl.InteriorPointOptions(
#         r_tol = 1.0e-8,
#         κ_init = 1.0e-6,
#         κ_tol = 2.0e-6,
#         diff_sol = true),
#     sim_opts = ContactControl.SimulatorOptions(warmstart = true))
#
# simulate!(sim, verbose = true)
#
# nz = num_var(model, env)
# nθ = num_data(model)
#
#
# ################################################################################
# # Mehrotra12
# # - linear       NO
# # - euclidean    NO
# # - conic        NO
# ################################################################################
#
# z2, θ2 = get_initialization(ref_traj, t)
# ip2 = mehrotra(z2, θ2,
# 	ix = linearization_var_index(model, env)[1],
# 	iy1 = linearization_var_index(model, env)[2],
# 	iy2 = linearization_var_index(model, env)[3],
# 	idyn = linearization_term_index(model, env)[1],
# 	irst = linearization_term_index(model, env)[2],
# 	ibil = linearization_term_index(model, env)[3],
#     idx_ineq = inequality_indices(model, env),
#     idx_soc = soc_indices(model, env),
# 	r! = s.res.r!,
#     rm! = s.res.rm!,
#     rz! = s.res.rz!,
#     rθ! = s.res.rθ!,
#     rz = s.rz,
#     rθ = s.rθ,
#     opts = Mehrotra12Options(
#         max_iter_inner=30,
#         r_tol=1e-8,
#         κ_tol=1e-8,
# 		# verbose=true
# 		))
# interior_point_solve!(ip2, z2, θ2)
# r2 = zeros(nz)
# ip2.methods.r!(r2, ip2.z, ip2.θ, 0.0)
# @testset "Mehrotra12 Nonlinear" begin
# 	@test norm(r2, Inf) < 1e-8
# 	@test abs(norm(r2, Inf) - 7.84e-9) < 1e-11
# 	@test ip2.iterations == 9
# end
