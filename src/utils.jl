
function scn(a::Number; digits::Int=1, exp_digits::Int=1)
	typeof(a) <: Float64 ? nothing : return nothing
end

function scn(a::Float64; digits::Int=1, exp_digits::Int=1)
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
	if m == 10.0
		m = 1.0
		e += 1
	end
    if digits == 0
        m = Int(floor(m))
		strm = string(m)
	else
		strm = string(m)
		is_neg = m < 0.
		strm = strm*"0"^max(0, 2+digits+is_neg-length(strm))
    end
    sgn = a >= 0 ? " " : ""
    sgne = e >= 0 ? "+" : "-"

	stre = string(abs(e))
	stre = "0"^max(0, exp_digits - length(stre)) * stre
    return "$sgn$(strm)e$sgne$(stre)"
end

function delta!(Δx::SizedArray{Tuple{nx},T,1,1}, x::SizedArray{Tuple{nx},T,1,1},
    x_ref::SizedArray{Tuple{nx},T,1,1}) where {nx,T}
    Δx .= x
    Δx .-= x_ref
    return nothing
end

function delta!(Δx::SizedArray{Tuple{nx},T,1,1}, x::AbstractArray{T},
    x_ref::AbstractArray{T}) where {nx,T}
    Δx .= x
    Δx .-= x_ref
    return nothing
end

function set!(s::SubArray, x::SizedArray{Tuple{nx},T,1,1}) where {nx,T}
    s .= x
    return nothing
end

function setminus!(s::SubArray, x::SizedArray{Tuple{nx},T,1,1}) where {nx,T}
    s .= -1.0 .* x
    return nothing
end


function convert_video_to_gif(video_file_path::String, output_path::String="output.gif";
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

function convert_frames_to_video_and_gif(filename, overwrite::Bool=true)
    MeshCat.convert_frames_to_video(
        homedir() * "/Downloads/$filename.tar",
        homedir() * "/Documents/video/$filename.mp4", overwrite=overwrite)

    convert_video_to_gif(
        homedir() * "/Documents/video/$filename.mp4",
        homedir() * "/Documents/video/$filename.gif", overwrite=overwrite)
    return nothing
end

function save_markdown(path::String, content::String; overwrite::Bool=true)
	mode = overwrite ? "w" : "a"
	io = open(path, mode)
	write(io, content)
	close(io)
	return nothing
end

function horizontal_line(n::Int)
	@assert n >= 1
	out = "|" * " --- |"^n
	return out
end

function content_line(c::AbstractVector)
	n = length(c)
	@assert n >= 1
	out = "|"
	for i = 1:n
		out *= string(c[i]) * "|"
	end
	return out
end

function save_expressions(expr::Dict{Symbol,Expr},
	path::AbstractString="expr.jld2"; overwrite::Bool=false)
	path = abspath(path)
    if isfile(path)
        if overwrite
            rm(path)
        end
        if !overwrite
            @warn "file exists -- not overwriting"
        end
    end
	@save path expr
	@info("Saved output as $path")
	return nothing
end

function load_expressions(path::AbstractString="expr.jld2")
	path = abspath(path)
	if !isfile(path)
		error("Could not find the input file $path")
	end
	@load path expr
	return expr
end

function module_dir()
	return joinpath(@__DIR__, "..")
end

function axes_pair_to_quaternion(n1, n2)
	if norm(n1 + n2, Inf) < 1e-5
		n2 = n2 + 1e-5ones(3)
	end

	reg(x) = 1e-20 * (x == 0) + x
	# provides the quaternion that rotates n1 into n2, assuming n1 and n2 are normalized
	n1 ./= reg(norm(n1))
	n2 ./= reg(norm(n2))
	n3 = skew(n1)*n2
	cθ = n1' * n2 # cosine
	sθ = norm(n3) # sine
	axis = n3 ./ reg(sθ)
	tanθhalf = sθ / reg(1 + cθ)
	q = [1; tanθhalf * axis]
	q /= norm(q)
	return Quaternion(q...)
end
