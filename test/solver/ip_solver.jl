using Random
using LinearAlgebra
# const ContactControl = Main

# include(joinpath(module_dir(), "src", "solver", "interior_point.jl"))
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
		ip = deepcopy(ip0)
	end
	e2 = @belapsed begin
		z, θ = get_initialization(ref_traj, t)
		ip = deepcopy(ip0)
		interior_point_solve!(ip, z, θ)
	end
	@show e1
	@show e2
	return e2 - e1
end

# function interior_point_profiling(ref_traj, t, ip0, N::Int = 100)
# 	@profiler begin
# 		for i = 1:N
# 			z, θ = get_initialization(ref_traj, t)
# 			ip = deepcopy(ip0)
# 			interior_point_solve!(ip, z, θ)
# 		end
# 	end
# 	return nothing
# end

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
	idx = OptimizationIndices(model, env),
	r! = s.res.r!,
    rz! = s.res.rz!,
    rθ! = s.res.rθ!,
    rz = s.rz,
    rθ = s.rθ,
    opts = InteriorPointOptions(
        max_iter=30,
        r_tol=1e-8,
        κ_tol=1e-8,
		# verbose=true
		))

interior_point_solve!(ip1)
r1 = zeros(nz)
s.res.r!(r1, ip1.z, ip1.θ, 0.0)
@testset "Interior Point Non Linear" begin
	@test norm(r1, Inf) < 1e-8
	@test abs(norm(r1, Inf) - 3.66e-9) < 1e-11
	@test ip1.iterations == 10
end
# e_ip = 1e6 * interior_point_timing(ref_traj, t, ip1)

################################################################################
# Test Interior Point on the linearized problem
################################################################################
im_traj1 = ImplicitTraj(ref_traj, s;
	κ = 1e-8,
	opts = InteriorPointOptions(
			κ_tol = 1.0e-8,
			r_tol = 1.0e-8,
			diff_sol = true,
			max_iter=100,
			ϵ_min=0.05,
			solver=:empty_solver,
			verbose=true
			))

z1, θ1 = get_initialization(ref_traj, t)
ip1 = deepcopy(im_traj1.ip[10])
interior_point_solve!(ip1, z1, θ1)
r1 = zeros(nz)
s.res.r!(r1, ip1.z, ip1.θ, 0.0)
@testset "Interior Point Linear" begin
	@test norm(r1, Inf) < 1e-8
end

# function convergence_plot()
# 	vio = [5e1, 4.4e-2, 5.2e-3, 8.2e-5, 3.9e-9]
# 	plt = plot(
# 	legend = false,
# 	xtickfontsize=18,
# 	ytickfontsize=18,
# 	xlabelfontsize=28,
# 	ylabelfontsize=28,
# 	xlabel = "iteration",
# 	ylabel = "violation",
# 	yaxis = :log,
# 	xlims = [0.8, length(vio) + 0.2],
# 	ylims = [6e-11, 1.5e2],
# 	)
# 	for i in eachindex(vio)
# 		plot!(plt, vio[1:i], color = :black, linewidth = 4.0)
# 		scatter!(plt, vio[1:i], color = :black, markersize = 7.0)
# 		display(plt)
# 		png(plt, joinpath(@__DIR__, "plot_$i.png"))
# 	end
# 	return
# end

# convergence_plot()
