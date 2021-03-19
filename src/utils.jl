# Run tests locally
using Pkg
Pkg.activate(@__DIR__)
# Pkg.activate(joinpath(@__DIR__, "test")

# Pkg.add("BenchmarkTools")
# Pkg.add("Colors")
# Pkg.add("FFMPEG")
# Pkg.add("ForwardDiff")
# Pkg.add("JLD2")
# Pkg.add("MeshCat")
# Pkg.add("Symbolics")
# Pkg.add("Parameters")
# Pkg.add("Plots")
# Pkg.add("Rotations")
# Pkg.add("StaticArrays")
# Pkg.add("LinearAlgebra")
# Pkg.add("Logging")
# Pkg.add("Random")
# Pkg.add("SparseArrays")
# Pkg.add("Test")
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#


function scn(a::Number; digits::Int=1)
	typeof(a) <: Float64 ? nothing : return nothing
end

function scn(a::Float64; digits::Int=1)
	@assert digits >= 0
    # a = m x 10^e
    if a == 0
        e = 0
        m = 0.0
    else
        e = Int(floor(log(abs(a))/log(10)))
        m = a*exp(-e*log(10))
    end

    m = round(m, digits=digits)
    if digits == 0
        m = Int(floor(m))
		strm = string(m)
	else
		strm = string(m)
		is_neg = m < 0.
		strm = strm*"0"^abs(2+digits+is_neg-length(strm))
    end
    sgn = a >= 0 ? " " : ""
    sgne = e >= 0 ? "+" : ""
    return "$sgn$(strm)e$sgne$e"
end

#
