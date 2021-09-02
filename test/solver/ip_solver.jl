using Random
using LinearAlgebra
# const ContactControl = Main

# include(joinpath(module_dir(), "src", "solver", "interior_point.jl"))
# include(joinpath(module_dir(), "src", "solver", "interior_point_base.jl"))
# include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))

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
# Test Interior Point on the full non linear problem
################################################################################

z1, θ1 = get_initialization(ref_traj, t)
ip1 = interior_point(z1, θ1,
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
    opts = InteriorPointOptions(
        max_iter_inner=30,
        r_tol=1e-8,
        κ_tol=1e-8,
		# verbose=true
		))
interior_point_solve!(ip1)
r1 = zeros(nz)
s.res.r!(r1, ip1.z, ip1.θ, 0.0)
@testset "Interior Point Non Linear" begin
	@test norm(r1, Inf) < 1e-8
	@test abs(norm(r1, Inf) - 2.06e-9) < 1e-11
	@test ip1.iterations == 8
end
# e_ip = 1e6 * interior_point_timing(ref_traj, t, ip1)

################################################################################
# Test Mehrotra on the full non linear problem
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
    opts = MehrotraOptions(
        max_iter_inner=30,
        r_tol=1e-8,
        κ_tol=1e-8,
		# verbose=true
		))
interior_point_solve!(ip2, z2, θ2)
r2 = zeros(nz)
ip2.methods.r!(r2, ip2.z, ip2.θ, 0.0)
@testset "Mehrotra Nonlinear" begin
	@test norm(r2, Inf) < 1e-8
	@test abs(norm(r2, Inf) - 7.84e-9) < 1e-11
	@test ip2.iterations == 9
end

# e_me = 1e6 * interior_point_timing(ref_traj, t, im_traj2.ip[t])

# @profiler for k = 1:3000
# 	z2, θ2 = get_initialization(ref_traj, t)
# 	interior_point_solve!(ip2, z2, θ2)
# end
# @elapsed for k = 1:3000
# 	z2, θ2 = get_initialization(ref_traj, t)
# 	interior_point_solve!(ip2, z2, θ2)
# end


################################################################################
# Test Interior Point on the linearized problem
################################################################################

im_traj1 = ImplicitTraj(ref_traj, s;
	κ = 1e-8,
	ip_type = :interior_point,
	opts = InteriorPointOptions(
			κ_tol = 1.0e-8,
			r_tol = 1.0e-8,
			diff_sol = true,
			max_iter_inner=100,
			ϵ_min=0.05,
			solver=:empty_solver,
			# verbose=true
			))

z1, θ1 = get_initialization(ref_traj, t)
ip1 = deepcopy(im_traj1.ip[10])
interior_point_solve!(ip1, z1, θ1)
r1 = zeros(nz)
s.res.r!(r1, ip1.z, ip1.θ, 0.0)
@testset "Interior Point Linear" begin
	@test norm(r1, Inf) < 1e-8
end

# e_ip = 1e6 * interior_point_timing(ref_traj, t, im_traj1.ip[t])

# @profiler for k = 1:5000
# 	z1, θ1 = get_initialization(ref_traj, t)
# 	ip1 = deepcopy(im_traj1.ip[10])
# 	interior_point_solve!(ip1, z1, θ1)
# end



################################################################################
# Test Mehrotra on the linearized problem
################################################################################

im_traj2 = ImplicitTraj(ref_traj, s;
	κ = 1e-8,
	ip_type = :mehrotra,
	opts = MehrotraOptions(
			κ_tol = 1.0e-8,
			r_tol = 1.0e-8,
			diff_sol = true,
			max_iter_inner=100,
			ϵ_min=0.05,
			solver=:empty_solver,
			# verbose=true
			))
z2, θ2 = get_initialization(ref_traj, t)
ip2 = deepcopy(im_traj2.ip[t])
interior_point_solve!(ip2, z2, θ2)
r2 = zeros(nz)
s.res.r!(r2, ip2.z, ip2.θ, 0.0)
@testset "Mehrotra Linear" begin
	@test norm(r2, Inf) < 1e-8
	@test (norm(r2, Inf) - 4.4e-9) < 1e-13
	@test ip2.iterations == 3
end





#
# lin = LinearizedStep(s, ref_traj.z[t], ref_traj.θ[t], ref_traj.κ[1])
# opts = MehrotraOptions(solver = :empty_solver)
# ip = mehrotra(
# 		 deepcopy(ref_traj.z[t]),
# 		 deepcopy(ref_traj.θ[t]),
# 		 idx_ineq = inequality_indices(model, env),
# 		 ix = linearization_var_index(model, env)[1],
# 		 iy1 = linearization_var_index(model, env)[2],
# 		 iy2 = linearization_var_index(model, env)[3],
# 		 idyn = linearization_term_index(model, env)[1],
# 		 irst = linearization_term_index(model, env)[2],
# 		 ibil = linearization_term_index(model, env)[3],
# 		 r! = r!,
# 		 rm! = rm!,
# 		 rz! = rz!,
# 		 rθ! = rθ!,
# 		 r  = RLin(s, lin.z, lin.θ, lin.r, lin.rz, lin.rθ),
# 		 rz = RZLin(s, lin.rz),
# 		 rθ = RθLin(s, lin.rθ),
# 		 v_pr = view(zeros(1,1), 1,1),
# 		 v_du = view(zeros(1,1), 1,1),
# 		 opts = opts)
#
# interior_point_solve!(ip)



# e_me = 1e6 * mehrotra_timing(ref_traj, t, im_traj2.ip[t])

# @profiler for k = 1:5000
# 	z2, θ2 = get_initialization(ref_traj, t)
# 	ip2 = deepcopy(im_traj2[10])
# 	interior_point_solve!(ip2, z2, θ2)
# end
# @elapsed for k = 1:3000
# 	z2, θ2 = get_initialization(ref_traj, t)
# 	interior_point_solve!(ip2, z2, θ2)
# end



#
# ################################################################################
# # Test Mehrotra on the linearized problem
# ################################################################################
# im_traj2 = ImplicitTraj(ref_traj, s;
# 	κ = 1e-8,
# 	ip_type = :mehrotra,
# 	opts = MehrotraOptions(
# 			κ_tol = 1.0e-8,
# 			r_tol = 1.0e-8,
# 			diff_sol = true,
# 			max_iter_inner=100,
# 			ϵ_min=0.05,
# 			solver=:empty_solver,
# 			# verbose=true
# 			))
#
# im_traj3 = ImplicitTraj(ref_traj, s;
# 	κ = 1e-8,
# 	ip_type = :interior_point,
# 	opts = InteriorPointOptions(
# 			κ_tol = 1.0e-8,
# 			r_tol = 1.0e-8,
# 			diff_sol = true,
# 			max_iter_inner=100,
# 			ϵ_min=0.05,
# 			solver=:empty_solver,
# 			# verbose=true
# 			))
#
# for t = 1:ref_traj.H
# 	z, θ = get_initialization(ref_traj, t)
# 	ip2 = deepcopy(im_traj2.ip[t])
# 	interior_point_solve!(ip2, deepcopy(z), deepcopy(θ))
# 	r2 = zeros(nz)
# 	s.res.r!(r2, ip2.z, ip2.θ, 0.0)
#
# 	ip3 = deepcopy(im_traj3.ip[t])
# 	interior_point_solve!(ip3, deepcopy(z), deepcopy(θ))
# 	r3 = zeros(nz)
# 	s.res.r!(r3, ip3.z, ip3.θ, 0.0)
# 	@test norm(r2, Inf) < 1e-8
# 	@test norm(r3, Inf) < 1e-8
# 	# @test norm((r2 - r3) ./ max.(r3, 1e-10), Inf) < 1e-1
# 	# @show scn.(ip2.z - ip3.z)
# 	# @show scn(norm(ip2.z - ip3.z, Inf))
# 	# @test norm((ip2.z - ip3.z) , Inf) < 1e-1
# 	@test ip2.iterations == ip3.iterations
# end
#
#
#
#
# im_traj2 = ImplicitTraj(ref_traj, s;
# 	κ = 1e-8,
# 	ip_type = :mehrotra,
# 	opts = MehrotraOptions(
# 			κ_tol = 2.0 * 1e-8,
# 			r_tol = 1.0e-8,
# 			diff_sol = true,
# 			max_iter_inner=100,
# 			ϵ_min=0.05,
# 			solver=:empty_solver,
# 			verbose=true
# 			))
# z2, θ2 = get_initialization(ref_traj, t)
# ip2 = deepcopy(im_traj2.ip[t])
# interior_point_solve!(ip2, z2, θ2)
# r2 = zeros(nz)
# s.res.r!(r2, ip2.z, ip2.θ, 0.0)
# @testset "Mehrotra Linear" begin
# 	@test norm(r2, Inf) < 1e-8
# 	@test (norm(r2, Inf) - 4.4e-9) < 1e-13
# 	@test ip2.iterations == 3
# end
#
# im_traj3 = ImplicitTraj(ref_traj, s;
# 	κ = 1e-8,
# 	ip_type = :interior_point,
# 	opts = InteriorPointOptions(
# 			κ_tol = 2.0 * 1e-8,
# 			r_tol = 1.0e-8,
# 			diff_sol = true,
# 			max_iter_inner=100,
# 			ϵ_min=0.05,
# 			max_ls = 1,
# 			solver=:empty_solver,
# 			verbose=true
# 			))
# z3, θ3 = get_initialization(ref_traj, t)
# ip3 = deepcopy(im_traj3.ip[t])
# interior_point_solve!(ip3, z3, θ3)
# r3 = zeros(nz)
# s.res.r!(r3, ip3.z, ip3.θ, 0.0)
# @testset "Mehrotra Linear" begin
# 	@test norm(r3, Inf) < 1e-8
# 	@test (norm(r3, Inf) - 4.4e-9) < 1e-13
# 	@test ip3.iterations == 3
# end
