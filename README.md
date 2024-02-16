# pTatinIO.jl
The aim of this package is to read files, produced with [pTatin3D](https://bitbucket.org/dmay/ptatin-total-dev.git), back into julia. 

In the future, we plan to expand this to be able to create initial model setups using the [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) package.

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

### Using it to load a timestep back into julia

We assume that you have used [ptatin3d-pyviztools](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) to process all the timesteps of a pTatin3D simulation, that your output data is in the directory`"test3"` relative to your current directory, and that you want to read the file `"step0-cellfields_compress.vts"`:
```julia
julia> using pTatinIO
julia> pTatin_CartData("step0-cellfields_compress.vts", Dir="test3/")
CartData 
    size    : (32, 32, 32)
    x       ϵ [ 3.90625 : 96.09375]
    y       ϵ [ 3.90625 : 96.09375]
    z       ϵ [ 3.90625 : 96.09375]
    fields  : (:region, :density, :viscosity, :e2, :RSR, :StrainRate, :T, :u)
  attributes: ["note"]
```
The output is in the [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) format.

### Reducing the filesize
One of the neat features of the `WriteVTK` package is that it writes compressed vtk files to disk, which can reduce the filesize by half.

```julia
julia> Compress_pTatin_Simulation("timeseries-cellfields.pvd", Dir="test3")
"timeseries-cellfields_compressed.pvd
```



