# ContactImplicitMPC.jl examples

This directory contains various ContactImplicitMPC.jl examples.
The `.jl` files in each subdirectory are meant to be processed using [Literate.jl](https://github.com/fredrikekre/Literate.jl).
During the documentation build process, the `.jl` files are converted to markdown
files that end up in the package documentation.

You can also generate Jupyter notebooks and run them locally by performing the following steps:

1. [install ContactImplicitMPC.jl](https://github.com/thowell/ContactImplicitMPC.jl)
2. [install IJulia](https://github.com/JuliaLang/IJulia.jl) (`add` it to the default project)
3. in the Julia REPL, run
   ```
   using Pkg
   Pkg.build("ContactImplicitMPC")
   using IJulia, ContactImplicitMPC
   notebook(dir=joinpath(dirname(pathof(ContactImplicitMPC)), "..", "examples"))
   ```

