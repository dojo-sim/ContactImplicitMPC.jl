# ContactImplicitMPC.jl examples

This directory contains examples of contact-implicit model-predictive control.
The `.jl` files in each subdirectory are meant to be processed using [Literate.jl](https://github.com/fredrikekre/Literate.jl).
During the build process, the `.jl` files are converted to notebooks. 

The initial installation of the package includes a number of pre-generated examples: [flamingo](examples/flamingo/flat.jl), [pushbot](examples/pushbot/push_recovery.jl), [hopper](examples/hopper/flat.jl), and [quadruped](examples/quadruped/flat.jl). You can generate Jupyter notebooks and run them locally by performing the following steps:

1. [install ContactImplicitMPC.jl](https://github.com/thowell/ContactImplicitMPC.jl)
2. [install IJulia](https://github.com/JuliaLang/IJulia.jl) (`add` it to the default project)
3. in the Julia REPL, run (do once)
   ```
   using Pkg
   Pkg.build("ContactImplicitMPC")
   ```
4. interact with notebooks
   ```
   using IJulia, ContactImplicitMPC
   notebook(dir=joinpath(dirname(pathof(ContactImplicitMPC)), "..", "examples"))
   ```

Additional examples can be run by first generating the models and simulations. Note that this may take 30-60 minutes.

5. in the Julia REPL, run (do once)
   ```
   include(joinpath(dirname(pathof(ContactImplicitMPC)), "..", "examples/generate_examples.jl"))
   ```