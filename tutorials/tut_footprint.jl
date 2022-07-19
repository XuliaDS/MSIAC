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
st = [true,false,true,false,true,false]
mesh_files = [create_2D_mesh(4,nX,nY,xBounds,yBounds,false,"m1.txt",structured=st[1],pert  = [0.2,0.2]),
              create_2D_mesh(4,nX,nY,xBounds,yBounds,false,"m2.txt",structured=st[2],pert = [0.1,0.1]),
              create_2D_mesh(3,nX,nY,xBounds,yBounds,false,"m3.txt",structured=st[3],pert  = [0.2,0.2]),
              create_2D_mesh(3,nX,nY,xBounds,yBounds,false,"m4.txt",structured=st[4],pert = [0.1,0.1]),
              create_2D_mesh(4,nX,nY,xBounds,yBounds.+[-1,1],true ,"m5.txt"),
              create_2D_mesh(3,nX,nY,xBounds,yBounds,true ,"m6.txt",structured=st[6])]

z = GZW("legendre",6)


for j= 2:2
    m    = "m"*string(j)*".txt"
    e    = shuffle([1:nX*nY;])[1]
    zI1  = shuffle([1:6;])[1]
    zI2  = shuffle([1:6;])[1]
    p    = shuffle([1,2,3])[1]
    zeta = [z.nodes[zI1],z.nodes[zI2]]#[z.nodes[zI1[i]],z.nodes[zI2[i]]]
    # we are going to mimmic a L2 field of degree p
    ff = create_analytic_field([p+1,p+1],["legendre","legendre"],m,"field.txt",type="ctt")
    # we now load our data for post-processing
    t2 = @elapsed data = load_data(m,ff, degree = p, structured=st[j]);  #Re-project data using a Pk basis
    f1 = test_filter(2*p+1,p+1,2*p+1,p+1,data,e,zeta,"line", track_footprint=true)
    f2 = test_filter(2*p+1,p+1,2*p+1,p+1,data,e,zeta,"tensor", track_footprint=true)
    @info " Convolution should be ~1 : Line = $f1 Tensor = $f2"
end
