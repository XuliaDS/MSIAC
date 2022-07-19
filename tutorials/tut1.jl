push!(LOAD_PATH, ".")
push!(LOAD_PATH, ENV["MSIAC"]) # tell Julia where to find MSIAC
using MSIAC   #(you can also write import MSIAC and then add module functions using MSIAC.function)
# we are going to plot a kernel of degree p = 2
   p = 2
   test_kernel(p)
   type_sol="sine"
# we are going to generate a triangular mesh
   nX = 10 ; nY = 10 ;
   xBounds = [0.0,2.0] ; yBounds = [0.0,2.0];
   periodic = true;
   mf = create_2D_mesh(3,nX,nY,xBounds,yBounds,periodic,"myTriMesh.txt", structured=true, pert = [0.0,0.0])

# we are going to mimmic a L2 field of degree p
   ff = create_analytic_field([p+1,p+1],["legendre","legendre"],mf,"mySineField.txt",type=type_sol)

# we now load our data for post-processing
   data = load_data(mf,ff, modalExp="Pk");  #Re-project data using a Pk basis

#checkout the mesh
        plot_mesh(data)
        #if you don't want to see the vertex & element map:
        plot_mesh(data, noIds = true)

# lets sample our data in a larger set of points
 update_phys!(data, "legendre",6,"legendre",6)
proj = get_field(data)
# lets look out our data
   plot_field(data, f=proj)

# we now filter our data
   @time filtered = filter_data(data, "line");

# lets look out our filtered data
   plot_field(data, f=filtered)

# compare the errors
 p1,p2 = get_evaluation_points(data)
 exact = get_exact_data(data, type_sol, p1,p2)
  @info " Linf error projection: $(maximum(abs.(exact .- proj))), filtered $(maximum(abs.(exact .- filtered)))"

rm(ff)
rm(mf)
