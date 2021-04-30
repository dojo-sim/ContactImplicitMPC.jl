function linearized_solver_timing_gain(model::ContactDynamicsModel)
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
	return gain, t_naive, t_efficient
end

model = get_model("hopper_2D")
model = get_model("hopper_3D")
model = get_model("flamingo")
model = get_model("pushbot")
model = get_model("quadruped")
linearized_solver_timing_gain(model)

# # Test

struct MD str end                                                      #hide
Base.show(io::IO, ::MIME"text/markdown", md::MD) = print(io, md.str)   #hide
function f()                                                           #hide
    io = IOBuffer()                                                    #hide
    for i in 1:2                                                       #hide
        println(io, "## test_$(i)")                                    #hide
    end                                                                #hide
    return MD(String(take!(io)))                                       #hide
end                                                                    #hide
md = f()

typeof(md.str)

save_markdown(joinpath(@__DIR__, "test.md"), "## end", overwrite=false)




function save_markdown(path::String, content::String; overwrite::Bool=true)
	mode = overwrite ? "w" : "a"
	io = open(path, mode)
	write(io, content)
	close(io)
	return nothing
end
