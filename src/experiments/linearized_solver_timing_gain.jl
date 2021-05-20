mutable struct TimingGain{T}
	name::String
	gain::T
	t_naive::T
	t_efficient::T
end

function linearized_solver_timing_gain(model::ContactModel)
	# Sizes
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b
	nx = nq
	ny = nb + 2nc
	nz = num_var(model)
	nθ = num_data(model)

	# Linearization point
	z0 = rand(nz)
	θ0 = rand(nθ)
	r0 = rand(nz)
	rz0 = zeros(nz,nz)
	rθ0 = zeros(nz,nθ)

	model.res.rz!(rz0, z0, θ0)
	model.res.rθ!(rθ0, z0, θ0)

	# Evaluation point
	z = rand(nz)
	θ = rand(nθ)
	κ = 1e-4
	κv = [κ]

	# Efficient version of r and rz
	r1 = RLin(model, z0, θ0, r0, rz0, rθ0)
	r!(r1, z, θ, κ)
	rz1 = RZLin(model, rz0)
	rz!(rz1, z)

	# Naive version of r and rz
	r2  = rand(nz)
	model.linearized.r!(r2, z, θ, κv, z0, θ0, r0, rz0, rθ0)
	rz2 = rand(nz, nz)
	model.linearized.rz!(rz2, z, rz0)

	# Test linear_solve!
	Δ = rand(nz)
	# benchmark doesn't update Δ
	t_efficient = @belapsed linear_solve!(Δ, rz1, r1)
	t_naive = @belapsed rz2 \ r2
	# Update Δ
	linear_solve!(Δ, rz1, r1)
	gain = t_naive / t_efficient
	@assert norm(Δ - rz2 \ r2, Inf) < 1e-10
	return TimingGain(string(model_name(model)), gain, t_naive, t_efficient)
end

function generate_markdown(tgs::Vector{<:TimingGain})
    io = IOBuffer()
	ncol = 4
	println(io, "# Linearized Solver Timing Gains")
	println(io, content_line(["model", "gain", "naive solver time (s)", "efficient solver time (s)"]))
	println(io, horizontal_line(ncol))
    for tg in tgs
		println(io, content_line([tg.name, scn(tg.gain), scn(tg.t_naive), scn(tg.t_efficient)]))
    end
	md = String(take!(io))
	return md
end

# Run experiment for models in names
# names = ["hopper_2D", "hopper_3D", "pushbot", "flamingo", "quadruped"]
# models = get_model.(names)
# tgs = linearized_solver_timing_gain.(models)
# content = generate_markdown(tgs)
# save_markdown(joinpath(@__DIR__, "linearized_solver_timing_gain.md"), content, overwrite=true)
