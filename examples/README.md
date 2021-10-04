# ContactImplicitMPC.jl examples

This directory contains examples of contact-implicit model-predictive control.
The `.jl` files in each subdirectory are meant to be processed using [Literate.jl](https://github.com/fredrikekre/Literate.jl).
During the build process, the `.jl` files are converted to notebooks and we generate models and simulation environments for each of the examples. Please note that the build time may take 30-60 minutes, but only needs to be performed once during the initial installation.

You can generate Jupyter notebooks and run them locally by performing the following steps:

1. [install ContactImplicitMPC.jl](https://github.com/thowell/ContactImplicitMPC.jl)
2. [install IJulia](https://github.com/JuliaLang/IJulia.jl) (`add` it to the default project)
3. in the Julia REPL, run
   ```
   using Pkg
   Pkg.build("ContactImplicitMPC")
   using IJulia, ContactImplicitMPC
   ```
4. interact with notebooks
   ```
   notebook(dir=joinpath(dirname(pathof(ContactImplicitMPC)), "..", "examples"))
   ```

