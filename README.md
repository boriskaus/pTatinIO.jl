# pTatinIO.jl
The aim of this package is to read files, produced with pTatin3D, back into julia (and at some stage create input files).

Note that it is very much work in progress at this stage.

### Installing
Install it through the julia package manager with:
```julia
julia>]
pkg> add https://github.com/boriskaus/pTatinIO.jl
```
Next, use the backspace to return to the julia `REPL` and type:
```julia
julia>using pTatinIO
```

### Using it

Assume that your output data is in the following directory, relative to your current directory
```julia
julia> step_dir="test3/step0/step0"
```
You can read the velocity data into julia with:
```julia
julia> using pTatinIO
julia> X,V = Read_Vel(step_dir)
```








