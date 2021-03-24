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
function convert_video_to_gif(video_file_path::AbstractString, output_path::AbstractString="output.gif";
    framerate::Int=30, start_time=0., duration=1e3, overwrite=false, width::Int=1080, height::Int=-2, hq_colors::Bool=false)
    output_path = abspath(output_path)

    if !isfile(video_file_path)
        error("Could not find the input file $video_file_path")
    end
    if isfile(output_path) && !overwrite
        error("The output path $output_path already exists. To overwrite that file, you can pass `overwrite=true` to this function")
    end

    mktempdir() do tmpdir
        # run(MeshCat.unpack_cmd(video_file_path, tmpdir, ".mp4", nothing)) # unpack the .tar file
        # cmd = ["-r", string(framerate), "-i", "%07d.png", "-vcodec", "libx264", "-preset", "slow", "-crf", "18"]
        color_map = hq_colors ?
            "[0:v] fps=$framerate, scale=$width:$height,split [a][b];[a] palettegen=stats_mode=single [p];[b][p] paletteuse=new=1" :
            "[0:v] fps=$framerate, scale=$width:$height,split [a][b];[a] palettegen [p];[b][p] paletteuse"
        cmd = ["-ss", string(start_time), "-t", string(duration), "-i", video_file_path, "-filter_complex", color_map]
        if overwrite
            push!(cmd, "-y")
        end
        push!(cmd, output_path)

        cd(tmpdir) do
            FFMPEG.exe(cmd...)
        end
    end
    @info("Saved output as $output_path")
    return output_path
end
