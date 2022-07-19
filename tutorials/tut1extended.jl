push!(LOAD_PATH, ".")
push!(LOAD_PATH, ENV["MSIAC"]) # tell Julia where to find MSIAC
#using Pkg
#Pkg.activate(normpath(ENV["MSIAC"]))
using MSIAC   #(you can also write import MSIAC and then add module functions using MSIAC.function)
using PyPlot
pygui(true)

@warn "Check that 2+ cores are used by typing: 1) on terminal: julia -p # {where #>1} OR 2) on Julia REPL: Threads.nthreads() {verify output is >1}."
@warn "If using Julia in VS Code, update number of threads in Settings>>Extensions>>Julia>> Julia: Num Threads. "

# Plot a kernel of degree p = 2
   p = 2
   @info "Kernel with degree $(p): (see plot)"
   test_kernel(p)
   #--- Consider a non-Symmetric Kernel either left or right shifted:
   #        shift {full left = -1.0, half left = - 0.5, full right = 1.0, half right = 0.5)}
   #        genSpline {1: with a general Spline, 0: without a general spline }
   @info "Non-symmetric kernel with degree $(p) shifted to the right with a general spline: (see plot)"
   test_kernel(p, shift = 1.0 , genSpline=1)
   @info "Non-symmetric kernel with degree $(p) shifted to the right without a general spline: (see plot)"
   test_kernel(p, shift = 1.0 , genSpline=0)

# Generate a Mesh:
   nX = 10 ; nY = 10 ;
   xBounds = [0.0,2.0] ; yBounds = [0.0,2.0];
   periodic = true;
   @info "~~~~~~~~ Different Mesh Files generated as named below: "
   #--- Uniform, Structured Quadrilateral mesh, type: 4 (first input), with mesh file name = "myQuadMesh.txt"
   quadFile = create_2D_mesh(4,nX,nY,xBounds,yBounds,periodic,"myQuadMesh.txt", structured=true, pert = [0.0,0.0])
   #--- Uniform, Structured Triangular mesh, type: 3, with mesh file name = "myTriMesh.txt"
   triFile = create_2D_mesh(3,nX,nY,xBounds,yBounds,periodic,"myTriMesh.txt", structured=true, pert = [0.0,0.0])
   #--- For Non-uniform, Non-structured Meshes indexed along vertical direction, consider other Options:
   #        structured {true or false},
   #        pert {perturb the uniform mesh by e.g. [deltaX, deltaY]=[0.4, 0.5]},
   #        yFirst {true=mesh indices run vertically, default: run horizontally}
   #.........Some Examples below, as named
   triFilePert = create_2D_mesh(3,nX,nY,xBounds,yBounds,periodic,"myTriMeshPert.txt", structured=true, pert = [0.4,0.5], yFirst=true)
   triFileUnstructPert = create_2D_mesh(3,nX,nY,xBounds,yBounds,periodic,"myTriMeshUnstructPert.txt", structured=false, pert = [0.4,0.5], yFirst=true)
   quadFilePert = create_2D_mesh(4,nX,nY,xBounds,yBounds,periodic,"myQuadMeshPert.txt", structured=true, pert = [0.2,0.4])
   quadFileUnstructPert = create_2D_mesh(4,nX,nY,xBounds,yBounds,periodic,"myQuadMeshUnstructPert.txt", structured=false, pert = [0.2,0.4])

# Visualize the Mesh Options:

   @info "$(nX)x$(nY) Triangular, Perturbed Mesh: (see plot)"
   mesh1 = load_mesh(triFilePert)
   plot_mesh(mesh1, noIds=true)

   @info "$(nX)x$(nY) Triangular, Perturbed, Unstructured Mesh: (see plot)"
   mesh2 = load_mesh(triFileUnstructPert)
   plot_mesh(mesh2, noIds=true)

   @info "$(nX)x$(nY) Quadrilateral, Perturbed Mesh: (see plot)"
   mesh3 = load_mesh(quadFilePert)
   plot_mesh(mesh3, noIds=true)

   @info "$(nX)x$(nY) Quadrilateral, Perturbed, Unstructured Mesh: (see plot)"
   mesh4 = load_mesh(quadFileUnstructPert)
   plot_mesh(mesh4, noIds=true)

# Select the Desired Mesh:
   meshFile = quadFile

# Generate an L2 field of degree p:
   type_sol = "sine"
   #--- Field with q=[p+1, p+1] data points along each direction of type Legendre
   #     with output file name = "mySineField.txt" for the Sine (default) analytic function
   fieldFile = create_analytic_field([p+1,p+1],["legendre","legendre"],meshFile,"mySineField.txt",type=type_sol)
   #--- Same Field with Lobatto data points
   fieldFile = create_analytic_field([p+1,p+1],["lobatto","lobatto"],meshFile,"mySineField.txt",type=type_sol)
   #--- For more Options:
   #        type {"sine" (sin(pi*x)cos(pi*x), "sinpi" (sin(pi*x), "poly_xy" (xy), ...; See eval_analytic_field in fields.jl for complete list)}
   #        time {for time-dependent analytic solutions, such as "burgers_sinXpY" (sin(pi*(x+y-2*time)))}
   type_sol = "sine"
   fieldFile = create_analytic_field([2*p+1,2*p+1],["legendre","legendre"],meshFile,"mySineField.txt",type=type_sol, time=1.0)

# Load Data for Post-Processing:
   data = load_data(meshFile,fieldFile);
   #--- For more Options:
   #        degree {select degree}
   #        modalExp {"Pk" for Pk basis, ...; default expansion is "legendre" for quads and "hierarchy" (modified legendre) for tris}
   #        structured {defaults as false}
   #        reverse {defaults as false}
   data = load_data(meshFile,fieldFile, modalExp="Pk", structured=false, reverse=false);
   data = load_data(meshFile,fieldFile, degree=p);

# Visualize the Mesh one More Time, with some more options:
   @info "Final Mesh File Selected with element ids: (see plot)"
   plot_mesh(data)
   #--- For more Options:
   #        noIds {defaults to false, set to true to remove the vertex & element map}
   #        blockPlot {defaults to true to show the plot}
   #        savePlot {filename to be saved as a .pdf image}
   @info "Final Mesh File Selected without element ids: (see plot)"
   plot_mesh(data, noIds=true, blockPlot=false)
   @info "Final Mesh File Selected without element ids saved to 'myMeshImg.pdf'"
   plot_mesh(data, noIds=true, blockPlot=false, savePlot="myMeshImg")

# Sample the Data in a larger set of points:
   #--- Sample the data field at 6 quadrature nodes  of the Legendre type for both the x and y directions
   update_phys!(data, "legendre",6,"legendre",6)
   #--- Get the field values
   proj = get_field(data) ;

# Visualize the Data:
   @info "Data field of type $(type_sol): (see plot)"
   plot_field(data, f=proj)
   #--- For more Options:
   #        showMesh {defaults to false, set to true to show the mesh}
   @info "Data field of type $(type_sol) with underlying mesh: (see plot)"
   plot_field(data, f=proj, showMesh=true)

# Filter the Data:
   #--- Filter the data with a Line Kernel, @time records the processing time
   @info "Processing time to filter data with Line filter: "
   @time filtered = filter_data(data, "line");
   #--- Filter the data with a Tensor Product Kernel
   @info "Processing time to filter data with Tensor Product filter: "
   @time filtered = filter_data(data, "tensor");
   #--- For more Options:
   #     parallel {defaults to true}
   #     rx {number of splines in the x direction}
   #     lx {order of splines in the x direction}
   #     ry {number of splines in the y direction}
   #     ly {order of splines in the y direction}
   #     dumpFile {default is "nofile", else: saves Fields to provided filename}
   #     theta {rotation in the kernel axis, defaults to "auto"}
   #     scaling {kernel scaling, if scaling=1 -> kernel scaling=1, if scaling !=1 -> kernel scaling=element size}
   #.....Comparing with and without multiple threads
   @info "Processing time to filter data with Line filter in parallel: "
   @time filtered = filter_data(data, "line", parallel=true);
   @info "Processing time to filter data with Line filter not in parallel: "
   @time filtered = filter_data(data, "line", parallel=false);
   #.....Considering different kernel options
   @info "Processing file saved under myDumpFile. Processing time to filter data with Line filter and designated parameters: "
   @time filtered = filter_data(data, "line", rx=2*p, lx=p, ry=2*p-1, ly=p, dumpFile="myDumpFile", theta="auto", scaling="auto");

   @info "Loading the filtered data from the generated file"
   data2 = load_data(meshFile,"myDumpFile.txt");
   filt2 = get_field(data2) ;
   @info " Checking solution is loaded correctly filt2 == filtered ?  $(filt2 == filtered) "

# Visualize the Filtered Data:
   @info "Filtered Data field of type $(type_sol) with underlying mesh: (see plot)"
   plot_field(data, f=filtered, showMesh=true)

# Compare Errors:
   #--- A bit of a cheating step here, where the same data is evaluated through eval_analytic_field to generate the "true" field
   #     Note: Above, this same step is built in the create_analytic_field used to generate fieldFile
   p1,p2 = get_evaluation_points(data)
   exact = get_exact_data(data, type_sol, p1,p2)
   size(proj)
   @info " Linf error projection: $(maximum(abs.(exact .- proj))), filtered $(maximum(abs.(exact .- filtered)))"

rm(fieldFile)
rm(meshFile)
