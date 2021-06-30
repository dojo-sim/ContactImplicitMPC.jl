function linearization_var_index(model::ContactModel, env::Environment)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)

	ny = 2nc + nb

	off = 0
	ix  = off .+ Vector(1:nq); off += nq
	iy1 = off .+ Vector(1:ny); off += ny
	iy2 = off .+ [Vector(nb .+ (1:nc)); Vector(1:nb); Vector(nb+nc .+ (1:nc))]; off += ny
	return ix, iy1, iy2
end

function linearization_term_index(model::ContactModel, env::Environment)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)

	ny = 2nc + nb
	# dyn = [dyn]
	# rst = [s1  - ..., ≡ ialt
	#        ... - η1 ,
	#        s2  - ...,]
	# bil = [γ1 .* s1 .- κ;
	#        b1 .* η1 .- κ;
	#        ψ1 .* s2 .- κ]
	idyn  = Vector(1:nq)
	irst1 = Vector(nq .+ (1:nc))
	irst2 = Vector(nq + nc .+ (1:nb))
	irst3 = Vector(nq + nc + nb .+ (1:nc))
	irst  = [irst1; irst2; irst3]
	ialt  = irst1
	ibil1 = Vector(nq + nb + 2nc .+ (1:nc))
	ibil2 = Vector(nq + nb + 3nc .+ (1:nb))
	ibil3 = Vector(nq + nb + 3nc + nb .+ (1:nc))
	ibil  = [ibil1; ibil2; ibil3]
	return idyn, irst, ibil, ialt
end

function get_bilinear_indices(model::ContactModel, env::Environment)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)

	terms = [SVector{nc,Int}(nq + nb + 2nc .+ (1:nc)),
			 SVector{nb,Int}(nq + nb + 3nc .+ (1:nb)),
			 SVector{nc,Int}(nq + nb + 3nc + nb .+ (1:nc))]

	vars = [[SVector{nc,Int}(nq .+ (1:nc)),           SVector{nc}(nq + 2nc + 2nb .+ (1:nc))], # γ1, s1
			[SVector{nb,Int}(nq + nc .+ (1:nb)),      SVector{nb}(nq + 2nc +  nb .+ (1:nb))], # b1, η
			[SVector{nc,Int}(nq + nc + nb .+ (1:nc)), SVector{nc}(nq + 3nc + 2nb .+ (1:nc))], # ψ, s2
			]
	return terms, vars
end

"""
    Structure holding the residual of the linearized residual r. This residual linearization is computed at
    a linearization point z0, θ0.
    The residual is decomposed into blocks to best exploit the structure of the problem.
         [rdyn ] -> Linearized dynamics constraints -> size = nx
    r  = |rrst | -> Rest of the linearized constraints -> size = ny
         [rbil ] -> Bilinear complementarity constraints -> size = ny
    The linear parts:
        [rdyn]
		|rrst| = [r0 + rz0*(z-z0) + rθ0*(θ-θ0)
		[....]
    The bilinear part:
        rbil = y1 .* y2 .- κ
"""
mutable struct RLin{T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}
    # Reference residual
    rdyn0::SVector{nx,T}
    rrst0::SVector{ny,T}
    rbil0::SVector{ny,T}
    # Residual
    rdyn::SVector{nx,T}
    rrst::SVector{ny,T}
    rbil::SVector{ny,T}

    # Reference residual jacobian rz0
    Dx::SMatrix{nx,nx,T,nxx}
    Dy1::SMatrix{nx,ny,T,nxy}
    Rx::SMatrix{ny,nx,T,nxy}
    Ry1::SMatrix{ny,ny,T,nyy}
    Ry2::SVector{ny,T}
    # Reference residual jacobian rθ0
    rθdyn::SMatrix{nx,nθ,T,nxθ}
    rθrst::SMatrix{ny,nθ,T,nyθ}
    rθbil::SMatrix{ny,nθ,T,nyθ}

    # Reference z0 and θ0
    x0::SVector{nx,T}
    y10::SVector{ny,T}
    y20::SVector{ny,T}
    θ0::SVector{nθ,T}
    # Values of z and θ
    x::SVector{nx,T}
    y1::SVector{ny,T}
    y2::SVector{ny,T}
    θ::SVector{nθ,T}

    # Indices
    nz::Int
    nθ::Int
    ix::SVector{nx,Int}
    iy1::SVector{ny,Int}
    iy2::SVector{ny,Int}
    iθ::SVector{nθ,Int}
    idyn::SVector{nx,Int}
    irst::SVector{ny,Int}
    ibil::SVector{ny,Int}
    ialt::SVector{nc,Int}

    # Altitude
    alt::SVector{nc,T}
    alt_zeros::SVector{nn,T}
end

function RLin(s::Simulation, z0::AbstractVector{T}, θ0::AbstractVector{T},
        r0::AbstractVector{T}, rz0::AbstractMatrix{T}, rθ0::AbstractMatrix{T}) where {T}

    model = s.model
    env = s.env

    nq = model.dim.q
    nc = model.dim.c
    nb = nc * friction_dim(env)
    nz = num_var(model, env)
    nθ = num_data(model)
    nx = nq
    ny = 2nc + nb
    nn = nc + nb

    # Terms
    off = 0
    idyn, irst, ibil, ialt = linearization_term_index(model, env)

    # Vars
    off = 0
    ix, iy1, iy2 = linearization_var_index(model, env)
    iθ = Vector(1:nθ)

    # Residual
    rdyn0 = SVector{nx,T}(r0[idyn])
    rrst0 = SVector{ny,T}(r0[irst])
    rbil0 = SVector{ny,T}(r0[ibil])
    rdyn = zeros(SVector{nx,T})
    rrst = zeros(SVector{ny,T})
    rbil = zeros(SVector{ny,T})

    # Matrices
    Dx = SMatrix{nx,nx,T,nx^2}(rz0[idyn,ix])
    Dy1 = SMatrix{nx,ny,T,nx*ny}(rz0[idyn,iy1])
    Rx = SMatrix{ny,nx,T,nx*ny}(rz0[irst,ix])
    Ry1 = SMatrix{ny,ny,T,ny^2}(rz0[irst,iy1])
    Ry2 = SVector{ny,T}(diag(rz0[irst,iy2]))
    y1 = SVector{ny,T}(diag(rz0[ibil,iy2]))
    y2 = SVector{ny,T}(diag(rz0[ibil,iy1]))
    rθdyn = SMatrix{nx,nθ,T,nx*nθ}(rθ0[idyn,iθ])
    rθrst = SMatrix{ny,nθ,T,ny*nθ}(rθ0[irst,iθ])
    rθbil = SMatrix{ny,nθ,T,ny*nθ}(rθ0[ibil,iθ])

    # Vectors
    x0  = SVector{nx,T}(z0[ix])
    y10 = SVector{ny,T}(z0[iy1])
    y20 = SVector{ny,T}(z0[iy2])
    θ0  = SVector{nθ,T}(θ0)
    x  = zeros(SVector{nx,T})
    y1 = zeros(SVector{ny,T})
    y2 = zeros(SVector{ny,T})
    θ  = zeros(SVector{nθ,T})

    return RLin{T,nx,ny,nθ,nx^2,nx*ny,ny*ny,nx*nθ,ny*nθ,nc,nn}(
        rdyn0,
        rrst0,
        rbil0,
        rdyn,
        rrst,
        rbil,

        Dx,
        Dy1,
        Rx,
        Ry1,
        Ry2,
        rθdyn,
        rθrst,
        rθbil,

        x0,
        y10,
        y20,
        θ0,
        x,
        y1,
        y2,
        θ,

        nz,
        nθ,
        SVector{nx,Int}(ix),
        SVector{ny,Int}(iy1),
        SVector{ny,Int}(iy2),
        SVector{nθ,Int}(iθ),
        SVector{nx,Int}(idyn),
        SVector{ny,Int}(irst),
        SVector{ny,Int}(ibil),
        SVector{nc,Int}(ialt),
        zeros(SVector{nc,T}),
        zeros(SVector{nc + nb,T}))
end

"""
    Structure holding the Jacobian of the linearized residual rz. This Jacobian linearization is computed at
    a linearization point z0, θ0.
    The Jacobian is decomposed into blocks to best exploit its sparsity and structure.
         [Dx  Dy1       0        ] -> Linearized dynamics constraints -> size = nx
    rz = |Rx  Ry1       Diag(Ry2)| -> Rest of the linearized constraints -> size = ny
         [0   Diag(y2)  Diag(y1) ] -> Bilinear complementarity constraints -> size = ny
          |      |         |
		  |		 |	       -> 1st set of variables associated with the bilinear constraints [γ1, b1, ψ] -> size = ny
          |      -> 2nd set of variables associated with the bilinear constraints [s1, η, s2] -> size = ny
		  -> Rest of the variables [q2] -> size = nx
    The constant parts are the ones that remains equals to the blocks in rz0:
        Dx, Dy1, Rx, Ry1, Ry2
    The variable parts:
        y1, y2.
	We define matrix A, it is used to initialize the interior point solver with
	a least square solution to the linearized constraints (i.e. excluding bilienar constraints).
	A = [Dx  Dy1       0        ] -> Linearized dynamics constraints -> size = nx
	    [Rx  Ry1       Diag(Ry2)] -> Rest of the linearized constraints -> size = ny
	δz = A' * inv(A * A') * δr
	We introduce the following notations:
	A' * inv(A * A') = [A1 A2]
					   |A3 A4|
	                   [A5 A6]
"""
mutable struct RZLin{T,nx,ny,nxx,nxy,nyy}
    # Reference residual jacobian rz0
    Dx::SMatrix{nx,nx,T,nxx}
    Dy1::SMatrix{nx,ny,T,nxy}
    Rx::SMatrix{ny,nx,T,nxy}
    Ry1::SMatrix{ny,ny,T,nyy}
    Ry2::SVector{ny,T}
    y1::SVector{ny,T}
    y2::SVector{ny,T}

    # Schur complement
    D::SMatrix{ny,ny,T,nyy}
    S::Schur{T,nx,ny,nxx,nxy,nyy}

	# Least Squares
	A1::SMatrix{nx,nx,T,nxx}
	A2::SMatrix{nx,ny,T,nxy}
	A3::SMatrix{ny,nx,T,nxy}
	A4::SMatrix{ny,ny,T,nyy}
	A5::SMatrix{ny,nx,T,nxy}
	A6::SMatrix{ny,ny,T,nyy}

    # Indices
    ix::SVector{nx,Int}
    iy1::SVector{ny,Int}
    iy2::SVector{ny,Int}
    idyn::SVector{nx,Int}
    irst::SVector{ny,Int}
    ibil::SVector{ny,Int}
end

function RZLin(s::Simulation, rz0::AbstractMatrix{T}) where {T}
    model = s.model
    env = s.env

    nq = model.dim.q
    nc = model.dim.c
    nb = nc * friction_dim(env)

    nz = num_var(model, env)
    nx = nq
    ny = 2nc + nb

    # Terms
    off = 0
    idyn, irst, ibil, ialt = linearization_term_index(model, env)

    # Vars
    off = 0
    ix, iy1, iy2 = linearization_var_index(model, env)

    # Fill the matrix blocks rz0s
    Dx = SMatrix{nx,nx,T,nx^2}(rz0[idyn,ix])
    Dy1 = SMatrix{nx,ny,T,nx*ny}(rz0[idyn,iy1])
    Rx = SMatrix{ny,nx,T,nx*ny}(rz0[irst,ix])
    Ry1 = SMatrix{ny,ny,T,ny^2}(rz0[irst,iy1])
    Ry2 = SVector{ny,T}(diag(rz0[irst,iy2]))
    y1 = SVector{ny,T}(diag(rz0[ibil,iy2]))
    y2 = SVector{ny,T}(diag(rz0[ibil,iy1]))

    # Schur complement
    D = Ry1 - Diagonal(Ry2 .* y2 ./ y1)
    M = [Dx Dy1;
         Rx D  ]
    S = Schur(M; n=nx, m=ny)

	# Least squares
	A = [Dx Dy1 zeros(nx, ny);
		 Rx Ry1 Diagonal(Ry2)]
	AA = A' * inv(A * A')
	A1 = SMatrix{nx,nx,T,nx*nx}(AA[1:nx,              1:nx])
	A2 = SMatrix{nx,ny,T,nx*ny}(AA[1:nx,              nx .+ (1:ny)])
	A3 = SMatrix{ny,nx,T,nx*ny}(AA[nx .+ (1:ny),      1:nx])
	A4 = SMatrix{ny,ny,T,ny*ny}(AA[nx .+ (1:ny),      nx .+ (1:ny)])
	A5 = SMatrix{ny,nx,T,nx*ny}(AA[nx + ny .+ (1:ny), 1:nx])
	A6 = SMatrix{ny,ny,T,ny*ny}(AA[nx + ny .+ (1:ny), nx .+ (1:ny)])

    return RZLin{T,nx,ny,nx^2,nx*ny,ny^2}(
        Dx,
        Dy1,
        Rx,
        Ry1,
        Ry2,
        y1,
        y2,
        D,
        S,
		A1,
		A2,
		A3,
		A4,
		A5,
		A6,
        SVector{nx,Int}(ix),
        SVector{ny,Int}(iy1),
        SVector{ny,Int}(iy2),
        SVector{nx,Int}(idyn),
        SVector{ny,Int}(irst),
        SVector{ny,Int}(ibil),
        )
end


"""
    Structure holding the Jacobian of the linearized residual rθ. This Jacobian linearization is computed at
    a linearization point z0, θ0.
"""
mutable struct RθLin{T,nx,ny,nθ,nxθ,nyθ}
    # Reference residual jacobian rθ0
	rθ0::Matrix{T}
	rθdyn0::SMatrix{nx,nθ,T,nxθ}
	rθrst0::SMatrix{ny,nθ,T,nyθ}
    rθbil0::SMatrix{ny,nθ,T,nyθ}

    # Indices
    idyn::SVector{nx,Int}
    irst::SVector{ny,Int}
    ibil::SVector{ny,Int}
end

function RθLin(s::Simulation, rθ0::AbstractMatrix{T}) where {T}
    model = s.model
    env = s.env

    nq = model.dim.q
    nc = model.dim.c
    nb = nc * friction_dim(env)

	nz = num_var(model, env)
    nθ = num_data(model)
    nx = nq
    ny = 2nc + nb

    # Terms
    off = 0

    idyn, irst, ibil, ialt = linearization_term_index(model, env)

    # Fill the matrix blocks rθ0s
	rθdyn0 = SMatrix{nx,nθ,T,nx*nθ}(rθ0[idyn,:])
	rθrst0 = SMatrix{ny,nθ,T,ny*nθ}(rθ0[irst,:])
	rθbil0 = SMatrix{ny,nθ,T,ny*nθ}(rθ0[ibil,:])

    return RθLin{T,nx,ny,nθ,nx*nθ,ny*nθ}(
        Matrix(rθ0),
		rθdyn0,
		rθrst0,
		rθbil0,
        SVector{nx,Int}(idyn),
        SVector{ny,Int}(irst),
        SVector{ny,Int}(ibil),
        )
end

"""
	Update the residual r.
"""
function r!(r::RLin{T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}, z::Vector{T}, θ::Vector{T}, κ::T,
        ) where {T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}
    r.x  = z[r.ix]
    r.y1 = z[r.iy1]
    r.y2 = z[r.iy2]
    r.θ  = θ[r.iθ]
    r.rdyn = r.rdyn0 + r.Dx*(r.x - r.x0) + r.Dy1*(r.y1 - r.y10)                         + r.rθdyn*(r.θ - r.θ0)
    r.rrst = r.rrst0 + r.Rx*(r.x - r.x0) + r.Ry1*(r.y1 - r.y10) + r.Ry2.*(r.y2 - r.y20) + r.rθrst*(r.θ - r.θ0) + SVector{ny}([r.alt; r.alt_zeros])
    r.rbil = r.y1 .* r.y2 .- κ
    return nothing
end

"""
	Update the residual rm.
"""
function rm!(r::RLin{T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}, z::Vector{T}, Δ::Vector{T}, θ::Vector{T}, κ::T,
        ) where {T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}
    r.x  = z[r.ix]
    r.y1 = z[r.iy1]
    r.y2 = z[r.iy2]
    r.θ  = θ[r.iθ]
    r.rdyn = r.rdyn0 + r.Dx*(r.x - r.x0) + r.Dy1*(r.y1 - r.y10)                         + r.rθdyn*(r.θ - r.θ0)
    r.rrst = r.rrst0 + r.Rx*(r.x - r.x0) + r.Ry1*(r.y1 - r.y10) + r.Ry2.*(r.y2 - r.y20) + r.rθrst*(r.θ - r.θ0) + SVector{ny}([r.alt; r.alt_zeros])
    r.rbil = r.y1 .* r.y2 + Δ[r.iy1] .* Δ[r.iy2] .- κ
    return nothing
end

function r_update!(r::RLin, r̄::RLin)
    r.x  = r̄.x
    r.y1 = r̄.y1
    r.y2 = r̄.y2
    r.θ  = r̄.θ
    r.rdyn = r̄.rdyn
    r.rrst = r̄.rrst
    r.rbil = r̄.rbil
    return nothing
end

"""
	Update the Jacobian rz, and update its Schur complement factorization.
"""
function rz!(rz::RZLin{T,nx,ny,nxx,nxy,nyy}, z::Vector{T}; reg::T = 0.0) where {T,nx,ny,nxx,nxy,nyy}
    # Unpack
    iy1 = rz.iy1
    iy2 = rz.iy2
    # id = rz.id
    Ry1 = rz.Ry1
    Ry2 = rz.Ry2

    # Update the matrix
    rz.y1 = z[iy1]
    rz.y2 = z[iy2]

	y1_reg = max.(rz.y1, reg)
	y2_reg = max.(rz.y2, reg)

    # update D in Schur complement
    rz.D = rz.Ry1 - Diagonal(rz.Ry2 .* y2_reg ./ y1_reg)
    # update Schur complement
    schur_factorize!(rz.S, rz.D)
    return nothing
end

"""
	Computes the search direction via a linear solve. The factorization has been done when updating rz,
	so this step should be very fast.
"""
function linear_solve!(Δ::Vector{T}, rz::RZLin{T,nx,ny,nxx,nxy,nyy},
        r::RLin{T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}; reg::T = 0.0) where {T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}
    # unpack
    rdyn = r.rdyn
    rrst = r.rrst
    rbil = r.rbil
    Ry2 = rz.Ry2
    y1 = rz.y1
    y2 = rz.y2

	y1_reg = max.(reg, y1)
	y2_reg = max.(reg, y2)

    u = rdyn
    v = rrst - Ry2 .*rbil ./ y1_reg
    schur_solve!(rz.S, u, v)
    Δ[rz.ix]  = rz.S.x
    Δ[rz.iy1] = rz.S.y
    Δ[rz.iy2] = (rbil .- y2_reg .* Δ[rz.iy1]) ./ y1_reg
    return nothing
end

"""
	Computes the Jacobian of the solution with respect to the parameters θ
	via a linear solve. The factorization has been done when updating rz,
	so this step should be ~ fast.
"""
function linear_solve!(δz::Matrix{T}, rz::RZLin{T,nx,ny,nxx,nxy,nyy},
	rθ::RθLin{T,nx,ny,nθ,nxθ,nyθ}; reg::T = 0.0) where {T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ}
    # unpack
	rθdyn = rθ.rθdyn0
	rθrst = rθ.rθrst0
	rθbil = rθ.rθbil0
	Ry2 = rz.Ry2
	y1  = rz.y1
	y2  = rz.y2
	ix  = rz.ix
	iy1 = rz.iy1
	iy2 = rz.iy2

	y1_reg = max.(reg, y1)
	y2_reg = max.(reg, y2)

	for i in eachindex(1:nθ)
		# We remark that rθbil is empty by construction of rθ.
		u = rθdyn[:,i]
		# v = rθrst[:,i] - Ry2 .* rθbil[:,i] ./ y1
		v = rθrst[:,i]
		schur_solve!(rz.S, u, v)
		@. δz[ix,i]  .= rz.S.x
		@. δz[iy1,i] .= rz.S.y
		δz[iy2,i] .= (rθbil[:,i] .- y2_reg .* δz[iy1,i]) ./ y1_reg
		# @. δz[iy2,i] .= .- y2 .* δz[iy1,i] ./ y1
	end
    return nothing
end

# Methods for the Interior Point solver
import LinearAlgebra.norm
function norm(r::RLin, t::Real)
	a = 0.0
	if t == Inf
		a = max(norm(r.rdyn, t),
				norm(r.rrst, t),
				norm(r.rbil, t))
	else
		a += norm(r.rdyn, t)
		a += norm(r.rrst, t)
		a += norm(r.rbil, t)
	end
	return a
end

function linear_solve!(solver::EmptySolver, δz::Matrix{T}, rz::RZLin{T,nx,ny,nxx,nxy,nyy},
	rθ::RθLin{T,nx,ny,nθ,nxθ,nyθ}; reg::T = 0.0) where {T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ}
	linear_solve!(δz, rz, rθ, reg = reg)
	return nothing
end

function linear_solve!(solver::EmptySolver, Δ::Vector{T}, rz::RZLin{T,nx,ny,nxx,nxy,nyy},
        r::RLin{T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}; reg::T = 0.0) where {T,nx,ny,nθ,nxx,nxy,nyy,nxθ,nyθ,nc,nn}
	linear_solve!(Δ, rz, r, reg = reg)
	return nothing
end


function update!(r::RLin{T}, z0::AbstractVector{T}, θ0::AbstractVector{T},
        r0::AbstractVector{T}, rz0::AbstractMatrix{T}, rθ0::AbstractMatrix{T}) where {T}
	idyn = r.idyn
	irst = r.irst
	ibil = r.ibil
	ix = r.ix
	iy1 = r.iy1
	iy2 = r.iy2

	# Reference residual
    r.rdyn0 = r0[idyn]
    r.rrst0 = r0[irst]
    r.rbil0 = r0[ibil]

    # Reference residual jacobian rz0
	r.Dx  = rz0[idyn, ix]
	r.Dy1 = rz0[idyn, iy1]
	r.Rx  = rz0[irst, ix]
	r.Ry1 = rz0[irst, iy1]
	r.Ry2 = diag(rz0[irst, iy2])

    # Reference residual jacobian rθ0
	r.rθdyn = rθ0[idyn,:]
	r.rθrst = rθ0[irst,:]
	r.rθbil = rθ0[ibil,:]

    # Reference z0 and θ0
    r.x0  = z0[ix]
    r.y10 = z0[iy1]
    r.y20 = z0[iy2]
    r.θ0  = θ0
	return nothing
end

function update!(rz::RZLin{T}, rz0::AbstractMatrix{T}) where {T}
	idyn = rz.idyn
	irst = rz.irst
	ibil = rz.ibil
	ix = rz.ix
	iy1 = rz.iy1
	iy2 = rz.iy2

    # Fill the matrix blocks rz0s
	rz.Dx  = rz0[idyn, ix]
	rz.Dy1 = rz0[idyn, iy1]
	rz.Rx  = rz0[irst, ix]
	rz.Ry1 = rz0[irst, iy1]
	rz.Ry2 = diag(rz0[irst, iy2])
	rz.y1  = diag(rz0[ibil, iy2])
	rz.y2  = diag(rz0[ibil, iy1])

	# Schur complement
	rz.D = rz.Ry1 - Diagonal(rz.Ry2 .* rz.y2 ./ rz.y1)
	M = [rz.Dx rz.Dy1;
	     rz.Rx rz.D  ]
 	nx = length(ix)
	ny = length(iy1)
	rz.S = Schur(M; n=nx, m=ny)
	return nothing
end

function update!(rθ::RθLin{T}, rθ0::AbstractMatrix{T}) where {T}
	# Fill the matrix blocks rθ0s
	rθ.rθdyn0 = rθ0[rθ.idyn,:]
	rθ.rθrst0 = rθ0[rθ.irst,:]
	rθ.rθbil0 = rθ0[rθ.ibil,:]
	return nothing
end

function bilinear_res(r::RLin, ibil)
    r.rbil
end

function least_squares!(z::Vector{T}, θ::AbstractVector{T}, r::RLin{T}, rz::RZLin{T}) where {T}
	# @warn "wrong"
	δθ = θ - r.θ0
	δrdyn = r.rdyn0 - r.rθdyn * δθ
	δrrst = r.rrst0 - r.rθrst * δθ

	δw1 = rz.A1 * δrdyn + rz.A2 * δrrst
	δw2 = rz.A3 * δrdyn + rz.A4 * δrrst
	δw3 = rz.A5 * δrdyn + rz.A6 * δrrst

	@. @inbounds z[r.ix]  .= r.x0  .+ δw1
	@. @inbounds z[r.iy1] .= r.y10 .+ δw2
	@. @inbounds z[r.iy2] .= r.y20 .+ δw3
	return nothing
end


# function r!(r::RLin{T}, z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
# 	r!(r, z, θ, κ)
# 	return nothing
# end

function rz!(rz::RZLin{T}, z::AbstractVector{T}, θ::AbstractVector{T}; reg::T = 0.0) where {T}
	rz!(rz, z, reg = reg)
	return nothing
end

function rθ!(rθ::RθLin{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where {T}
	return nothing
end
