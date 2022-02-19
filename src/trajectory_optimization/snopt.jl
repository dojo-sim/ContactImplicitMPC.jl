"""
    set path to SNOPT MathOptInterface interface.
        suggested: https://github.com/thowell/SNOPT7.jl
            follow installation instructions from README
        or: modify official repository (https://github.com/snopt/SNOPT7.jl)
            to have increased workspace size.

    SNOPT license is required.
"""
function include_snopt(path = joinpath(dirname(pwd()),"SNOPT7.jl/src/SNOPT7.jl"))
    include(path)
end
