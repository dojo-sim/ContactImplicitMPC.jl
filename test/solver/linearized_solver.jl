@testset "Linearized Solver" begin
	# model = ContactControl.get_model("quadruped")
	model = ContactControl.get_model("hopper_2D")

	# Sizes
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b
	nx = nq
	ny = nb + 2nc
	nz = ContactControl.num_var(model)
	nθ = ContactControl.num_data(model)

	# Indices
	rz_ = ContactControl.RZLin(model, rand(nz,nz))
	ix = rz_.ix
	iy1 = rz_.iy1
	iy2 = rz_.iy2
	idyn = rz_.idyn
	irst = rz_.irst
	ibil = rz_.ibil

	z0 = rand(nz)
	θ0 = rand(nθ)
	r0 = rand(nz)
	rz0 = zeros(nz,nz)
	rθ0 = zeros(nz,nθ)
	model.res.rz!(rz0, z0, θ0,  ContactControl.NoCache())
	model.res.rθ!(rθ0, z0, θ0,  ContactControl.NoCache())

	z = rand(nz)
	θ = rand(nθ)
	κ = 1e-4
	κv = [κ]
	rz1 = ContactControl.RZLin(model, rz0)
	r1 = ContactControl.RLin(model, z0, θ0, r0, rz0, rθ0)

	# Test rz!
	rz!(rz1, z)
	rz2 = rand(nz, nz)
	model.linearized.rz!(rz2, z, rz0)

	@test norm(rz1.Dx  - rz2[idyn, ix],  Inf) < 1e-10
	@test norm(rz1.Dy1 - rz2[idyn, iy1], Inf) < 1e-10
	@test norm(rz1.Rx  - rz2[irst, ix],  Inf) < 1e-10
	@test norm(rz1.Ry1 - rz2[irst, iy1],  Inf) < 1e-10
	@test norm(rz1.Ry2 - diag(rz2[irst, iy2]),  Inf) < 1e-10
	@test norm(rz1.y2  - diag(rz2[ibil, iy1]),  Inf) < 1e-10
	@test norm(rz1.y1  - diag(rz2[ibil, iy2]),  Inf) < 1e-10

	# Test r!
	r1 = ContactControl.RLin(model, z0, θ0, r0, rz0, rθ0)
	ContactControl.r!(r1, z, θ, κ)
	r2  = rand(nz)
	model.linearized.r!(r2, z, θ, κv, z0, θ0, r0, rz0, rθ0)
	@test norm(r2[r1.idyn] - r1.rdyn, Inf) < 1e-10
	@test norm(r2[r1.irst] - r1.rrst, Inf) < 1e-10
	@test norm(r2[r1.ibil] - r1.rbil, Inf) < 1e-10

	# Test linear_solve!
	Δ = rand(nz)
	ContactControl.linear_solve!(Δ, rz1, r1)
	@test norm(Δ - rz2\r2, Inf) < 1e-10

	# @benchmark ContactControl.rz!(rz1, z)
	# @benchmark ContactControl.r!(r1, z, θ, κ)
	# @benchmark ContactControl.linear_solve!(Δ, rz1, r1)

	# Test linear_solve! for matrices
	model.res.rθ!(rθ0, z, θ,  NoCache())
	@test norm(rθ0[ibil,:], Inf) < 1e-10
	# plot(Gray.(1e10*abs.(rθ_[idyn,:])))
	# plot(Gray.(1e10*abs.(rθ_[irst,:])))
	# plot(Gray.(1e10*abs.(rθ_[ibil,:])))

	rz1 = ContactControl.RZLin(model, rz0)
	rθ1 = ContactControl.RθLin(model, rθ0)
	δz1 = rand(nz,nθ)
	ContactControl.linear_solve!(δz1, rz1, rθ1)
	@test norm(δz1 - (rz0\rθ0), Inf) < 1e-10
	# @benchmark ContactControl.linear_solve!(δz1, rz1, rθ1)
end
