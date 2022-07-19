push!(LOAD_PATH, ".")
push!(LOAD_PATH, ENV["MSIAC"]) # tell Julia where to find MSIAC

using MSIAC   #(you can also write import MSIAC and then add module functions using MSIAC.function)
using PyPlot, PyCall, Random

# we are going to plot a kernel of degree p = 2
# we are going to generate a triangular mesh
periodic =true;
st  = "sine"
mIO = "mymesh.txt"
fIO = "mySine.txt"
nX = 12; nY = 12 ;
xBounds = [0.0,2.0] ; yBounds = [0.0,2.0];

mf = create_2D_mesh(4,nX,nY,xBounds,yBounds,false,"m1.txt",structured=true)
z = GZW("legendre",6)

e = 66

for p = 3:3
    zeta = [0.0,0.0]
    # we are going to mimmic a L2 field of degree p
    ff = create_analytic_field([p+1,p+1],["legendre","legendre"],mf,"field.txt",type="ctt")
    # we now load our data for post-processing
    t2 = @elapsed data = load_data(mf,ff, degree = p, structured=true);  #Re-project data using a Pk basis
    #plot_mesh(data)
    f1 = test_filter(2*p+1,p+1,2*p+1,p+1,data,e,zeta,"line", track_footprint=true)
    f2 = test_filter(2*p+1,p+1,2*p+1,p+1,data,e,zeta,"tensor", track_footprint=true)
    @info " Convolution should be ~1 : Line = $f1 Tensor = $f2"

end
