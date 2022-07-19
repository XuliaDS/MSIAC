
1) Delete all dependencies from your Project.toml file (everything under the line [deps])

Launch julia
julia --project=~/path_to_MSIAC 

--OR--

julia
julia> press ]
(pkg)> activate .

julia > include("MSIAC_pkg.jl") #this script simply adds manually to your local Julia all the packages that MSIAC depends on 

press ]

instantiate

build

***NOTE****: WriteVTK and VTKDataIO packages show some warnings during pre-compilation. However, they still work and allow reading .vtu & .vtk files correctly. 

press backspace

---- exit Julia

export JULIA_NUM_THREADS=8 (check your actual number of cores)


