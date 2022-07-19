#=
*
*      MSIAC: a Julia module designed to post-process data with SIAC filters
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
=#

export update_phys!, update_modes!, update_fields!, load_data,load_field, create_analytic_field, load_modes, load_modes_file, get_exact_data, eval_analytic_field,
       update_points!, get_evaluation_points, get_evaluation_weights, get_field

using  WriteVTK, VTKDataIO, PyCall, PyPlot
export plot_field,write_vtk, write_vtk_centres

get_exact_data(mesh::MESH, type::String, zw1::vecD, zw2::vecD; time = 0.0::dbl) = reshape(eval_analytic_field(mesh_local_to_global(mesh,zw1,zw2),type, time = time), length(zw1) * length(zw2), mesh.N)
get_exact_data(data::DATAFIELD, type::String, zw1::vecD, zw2::vecD; time = 0.0::dbl) = get_exact_data(data.mesh, type, zw1,zw2,time = time)
function eval_analytic_field(t::Union{matD,vecD},type::String;time=0.0::dbl)
	if type == "sine"
		return sinpi.(t[1,:]) .* cospi.(t[2,:])
	elseif type == "sinex"
		return sinpi.(t[1,:])
	elseif type == "ctt"
		return fill(1.0, length(t[1,:]))
	elseif type == "poly_xy"
		return t[2,:].* t[1,:]
	elseif type == "poly_x+y"
		return t[2,:].+ t[1,:]
	elseif type == "poly_xy2"
		return t[1,:] .* t[2,:].^2
	elseif type == "poly_x2y2"
		return t[1,:].^2 .* t[2,:].^2
	elseif type == "poly_xy3"
		return t[1,:] .* t[2,:].^3
	elseif type == "burgers_sinXpY"
		return sinpi.(t[1,:] .+ t[2,:] .- 2.0 * time)
	elseif type == "exp"
		return exp.(.-t[1,:].^2 .-t[2,:].^2)
	else
		return sinpi.(t[1,:].+ t[2,:])
	end
end

"""
Given a field file, reconstructs a field solution in modal form stored in a FIELD structure

Parameters:
===========
mesh file
field file

Optional:
=========
modalExp:  expansion type (default: legendre for quads and hierarchy for tris (modified legendre))
degree: polynomial degree for the expansion
reverse: by default (false), y runs first and x second. Reverse = true to switch it
Returns: a FIELD data structure
"""
function load_field(fileF::String; degree = -1, modalExp = "legendre",reverse= false)
    ## START LOADING MESH INFO
    #############################
	ext = split(fileF, ".")[end]
	if ext == "vtu" || ext == "vtk"
		return load_vtu(fileF)
	end
    data = readdlm(fileF)
	F, fInfo, fType = find_data("\$Fields", data)
	(F == nothing) && (error(" No field data !! ") && return nothing)
	n    = F[1]
	Q    = F[2]
	type = string(F[3])
	zwS1 = GZW(type, Q)
	Q    = F[4]
	type = string(F[5])
	zwS2 = GZW(type, Q)
	p    = degree == -1 ? Q - 1 : degree
	 # ASSUMINE ONLY ONE FIELD (SOLUTION) IS GIVEN
	bt   = fType[1] == 3 && modalExp == "legendre" ? "hierarchy" : modalExp
	return FIELD(bt,p, fInfo[1] == 3, zw = [zwS1,zwS2], f = hcat(fInfo...),reverse=reverse)
end

"""
Given a mesh & field file (ONE file if .vtk .vtu), creates a DATAFIELD structure storing the mesh and reconstructing a modal approximation from field points

Parameters:
===========
mesh file
field file
** one single .vtu (.vtk) file

Optional:
=========
structured: if mesh is structured, indicate true (faster computation)
modalExp:  expansion type (default: legendre for quads and hierarchy for tris (modified legendre))
degree: polynomial degree for the expansion
reverse: by default (false), y runs first and x second. Reverse = true to switch it
Returns: a DATAFIELD  structure
"""
function load_data(f1::String, f2=nothing::Union{String,Nothing}; structured=false,
	               degree = -1, modalExp = "legendre",reverse= false, structrued = false)
	ext = split(f1, ".")[end]
	if ext == "vtu" || ext == "vtk"
		return load_vtu(f1)
	end
	mesh  = load_mesh(f1; structured = structured)
	field = load_field(f2, degree = degree, modalExp = modalExp, reverse=reverse )
	return DATAFIELD([field], mesh)
end

"""
Creates a field file based on an analytic solution (for testing)
Parameters:
===========
q: number of data points in each direction[qx,qy]
eQ: type of data points in each direction (legendre, lobatto, etc)
mesh: the mesh file
name: name of the output file

Optional:
=========
type: analytic function (see eval_analytic_field) default: sine
time: for time dependent analytic solutions
"""
function create_analytic_field(q::Vector{Int}, eQT::Vector{String}, mesh::Union{MESH,String}, f_name::String;
	                            type = "sine"::String, time = 0.0::dbl)
	zw1 = GZW(eQT[1],q[1])
	zw2 = GZW(eQT[2],q[2])
	if typeof(mesh) == String
		mesh = load_mesh(mesh)
	end
    open(f_name, "w") do io
		println(io,"\$Fields")
		println(io,mesh.N," ",q[1]," ",eQT[1]," ",q[2]," ",eQT[2])
		for j = 1:mesh.N
			j_nodes = mesh.elmts[j].nodes
			print(io,j," ",length(j_nodes)," ")
			#t       = local_to_global(zw1.nodes,zw2.nodes, hcat(getfield.(m_nodes[j_nodes],:xyz)...))
			verts   = hcat(getfield.(mesh.nodes[j_nodes],:xyz)...)
			t,dumb  = standard_region(zw1.nodes,zw2.nodes, verts)#hcat(getfield.(m_nodes[j_nodes],:xyz)...))
			tt      = zeros(size(t))
			tt[1,:] = t[1,:]
			tt[2,:] = t[2,:]
			f       = eval_analytic_field(tt,type,time=time)
			[ print(io,dbl(f[j])," ") for j = 1:length(f)]
			print(io,"\n")
		end
		println(io,"\$EndFields")
		end
	return f_name
end

load_modes(q::Int, modes::matD, btype::String, collapse::Bool,name=""::String) = FIELD(btype,q-1, collapse,modes, name)

"""
Given a modes file, creates a FIELD data structure

Parameters:
===========
p = degree
m = total modes (x-dir + y-dir)
n = elements
btype = basis type corresponding to modes
collapse = true if triangular mesh
f = file name
Optional:
=========
name: name the field variable (pressure, u, etc)
Returns: a FIELD data structure
"""
function load_modes_file(p::Int, m::Int, n::Int, btype::String, collapse::Bool,f::String;name=""::String)
	data  = readdlm(f)
	modes = dbl.(reshape(data, m , n))
	 # ASSUMINE ONLY ONE FIELD (SOLUTION) IS GIVEN
	return FIELD(btype,p,collapse,m = modes)
end

function load_vtu(fname::String)
    ext = split(fname, ".")[end]
    if ext != "vtu" && ext != "vtk"
        @error " Invalid data file  $fname needs extension .vtu or .vtk"
        return nothing
    end
    gridreader = read_vtk(fname)
	coords = gridreader.point_coords
	conne  = gridreader.cell_connectivity
	# Elements map: adjacent nodes + coordinates
	vDict   = Dict{Int, Vector{Int}}()
	eDict   = Dict{Vector{Int},Vector{Int}}()
	nE      = length(conne)
	nV      = length(coords[1,:])
	nm      = [1:nV;]
	cXY     = zeros(3) ; cX = zeros(3) ; cY = zeros(3)
	m_nodes = Array{NODE, 1}(undef,nV)
	pdata   = gridreader.point_data
	fnames  = collect(keys(pdata))
	nf      = length(fnames)
	zw      = GZW("Lobatto",2)
	p       = 1
	nv      = length(conne[1])
	fields  = [FIELD("hierarchy",p,nv == 3,zw=[zw,zw],n = fnames[j],f = zeros(nv,nE), m = zeros((p+1)^2,nE), reverse=true) for j = 1:nf]
	for j = 1:nE
		ele = conne[j]
		ca  = length(ele) == 4 ? [1,2,4,3,1] : [1,2,3,1]
		for f = 1 : nf
			phi = pdata[fnames[f]]
			fields[f].f[:,j] = phi[ele]  # x runs first
		end
		for f = 1:length(ele)
			vp = [ele[ca[f]], ele[ca[f+1]]]
			vDict = add_key(vDict, vp[1],j)
			eDict = add_key(eDict, sort(vp), j)
			nm[vp[1]] == 0 && continue
			m_nodes[vp[1]] = NODE(coords[:,vp[1]], unique!(vDict[vp[1]]), [vp[1]])
			cj = coords[:,vp[1]]
			if j == 1
				cXY .= cj ; cY .= cj ; cX .= cj
			else
				all(cj.< cXY) && (cXY.= cj)
				cj[1] > cX[1] && (cX[1] = cj[1])
				cj[2] > cY[2] && (cY[2] = cj[2])
			end
			nm[vp[1]] = 0
		end
		[update_element_modes!(fields[f], [zw,zw], j) for f = 1:nf]
	end
	m_elmts = Array{ELEMENT, 1}(undef, nE)
	for j = 1:length(conne)
		nj = conne[j]
		le = length(nj)
		cE  = le == 4 ? [1,2,4,3,1] : [1,2,3,1]
		vj  = [m_nodes[nj[cE[k+1]]].xyz .- m_nodes[nj[cE[k]]].xyz for k = 1:le]
		hj  = norm.(vj)
		area = norm(cross(vj[1], vj[2]))
		le == 3 && (area *= 0.5)
		m_elmts[j] = ELEMENT(nj[cE[1:le]], zeros(Int,(le,2)), le,hj, area)
		for k = 1 : le
			nj1 = sort(nj[cE[k:k+1]])
			idx = findfirst(x->x != j, eDict[nj1])
			m_elmts[j].neigh[k,1] = idx != nothing ? eDict[nj1][idx] : 0
			m_elmts[j].neigh[k,2] = idx != nothing ? 1 : 0 ## -1 would have been periodic
		end
	end
	axis  = hcat([(cX .- cXY) ./ norm( cX .- cXY),
		  		  (cY .- cXY) ./ norm( cY .- cXY)]...)'

	mesh = MESH(m_nodes,m_elmts,[[0]], false, false, axis, nE, fname)
	return DATAFIELD(fields,mesh)
end

function eval_basis(basis::EXPANSION, z1::vecD, z2::vecD)
	P     = zeros(length(z1),length(basis.m1))
	if basis.type == "pk" || basis.type == "Pk" || basis.type == "PK"
		for j = 1:length(basis.m1)
	    	c       = (0.5)^(basis.m1[j] + 1 / 2) * sqrt((2 * basis.m1[j] + 1) * (basis.m1[j] + basis.m2[j] + 1))
	        JYvals  = poly_jacobi(basis.m2[j],2 * basis.m1[j] + 1, 0.0)
			JXvals  = poly_legendre(basis.m1[j])
			P[:, j] = c .* JXvals.(z1) .* JYvals.(z2) .* ((1 .- z2).^(basis.m1[j]))
	    end
	else
		[P[:,j] = basis.J[basis.m1[j]].(z1) .* basis.J[basis.m2[j]].(z2) for j = 1:length(basis.m1)]
	end
	return P
end

export eval_basis

"""
Updates the field(s) in data evaluated at the quadrature points (zw)
Parameters:
===========
data: a DATAFIELD / FIELD structure
zw: a 2D vector of GZW points

Optional:
=========
idx: if "all" updates all the fields in data, otherwise specify a particular field (p,u,v,etc)
"""
function update_phys!(data::Union{DATAFIELD,FIELD},q1::String, n1::Int,q2::String,n2::Int;idx = "all"::String)
	if typeof(data) == FIELD
		data.zw[1] = GZW(q1,n1); data.zw[2] = GZW(q2,n2)
		data.f = phys_values(data.basis, data.modes, data.zw, data.reverse)
	else
		piv = idx != "all" ? [idx] : getfield.(data.fields,:name)
		for j in piv
			k = data.fdict[j]
			data.fields[k].zw[1] = GZW(q1,n1); data.fields[k].zw[2] = GZW(q2,n2)
			data.fields[k].f = phys_values(data.fields[k].basis, data.fields[k].modes, data.fields[k].zw, data.fields[k].reverse)
		end
	end
end


"""
Given a field name (u,v,w etc) returns the data stored for that variable
"""
get_field(data::DATAFIELD, idx=""::String) = data.fields[data.fdict[idx]].f

function phys_at_element(data::FIELD,z1::vecD,z2::vecD,ele::Int)
	phi = eval_basis(data.basis, z1, z2)
	return [sum(data.modes[:,ele] .* phi[j,:]) for j = 1 : length(z1)*length(z2)]
end

function update_element_modes!(field::FIELD, zw::Vector{GZW}, e::Int)
	# set mass matrix
	f = field.f[:,e]
	Q = length(f)
	i1,i2 = tensor_index(zw[1].Q,zw[2].Q, field.reverse)
	if (Q != length(i1))
		@error " Inside modes_f :: MISMATCH BETWEEN WEIGHTS AND FIELD $Q != $(zw[1].Q * zw[2].Q)"
		return nothing
	end
	# compute (f, phi_j) for j = 1 : order
	w12  = zw[1].weights[i1] .* zw[2].weights[i2]
	phi  = eval_basis(field.basis,zw[1].nodes[i1],zw[2].nodes[i2])
	if field.basis.type == "Pk" || field.basis.type == "pk" || field.basis.type == "Pk"
		field.modes[:,e] = [sum(f .* w12 .* phi[:, j] .* 0.5 .* (1 .- zw[2].nodes[i2])) for j = 1:length(data.m1)]
	else
		field.modes[:,e] = field.basis.mass * [sum(f .* phi[:,j] .* w12) for j = 1 : length(phi[1,:])]
	end
end

"""
Updates the mode in data from the field values (also in data)
"""
update_modes!(data::FIELD) = (data.modes = modal_values(data.basis, data.f, data.zw, data.reverse))


"""
Physical (field) values given an element an its local coordinates (in reference system)
Parameters:
===========
DATAFIELD structure or FIELD structure
e: the element ID
z: points coordinates in reference system (chi: [-1,1]^2)

Returns:
========
Field values at point
"""
function evaluate_point(data::Union{DATAFIELD,FIELD},e::Int,z::vecD, idx="all"::String)
	if typeof(data) == FIELD
		modes  = data.modes[:,e]
		P      = eval_basis(data.basis,[z[1]],[z[2]])
		return sum(modes .* P[1,:])
	else
		piv = idx != "all" ? [idx] : getfield.(data.fields,:name)
		val = zeros(length(piv))
		for j = 1 : length(piv)
			k      = data.fdict[piv[j]]
			modes  = data.fields[k].modes[:,e]
			P      = eval_basis(data.basis,[z[1]],[z[2]])
			val[j] = sum(modes .* P[1,:])
		end
		return val
	end
end


"""
Returns the evaluation points (nodes) stored in data
Parameters:
===========
a FIELD structure
"""
get_evaluation_points(data::FIELD) = data.zw[1].nodes, data.zw[2].nodes
get_evaluation_points(data::DATAFIELD, idx=""::String) = get_evaluation_points(data.fields[data.fdict[idx]])

"""
Returns the evaluation weights stored in data
Parameters:
===========
a FIELD structure
"""
get_evaluation_weights(data::FIELD) = data.zw[1].weights, data.zw[2].weights
get_evaluation_weights(data::DATAFIELD, idx=""::String) =  get_evaluation_weights(data.fields[data.fdict[idx]])


plot_field(data::DATAFIELD;idx=""::String,f=matD(undef,0,0)::matD, showMesh=false::Bool) = plot_field(data.mesh, data.fields[data.fdict[idx]],f=f,showMesh = showMesh)

	function plot_field(msh::MESH,fld::FIELD;f=matD(undef,0,0)::matD, showMesh=false::Bool)
		ff = isempty(f) ? fld.f : f
		fig  = plt.figure(1)
		ax   = fig.add_subplot(111,projection="3d")

		t    = mesh_local_to_global(msh,fld.zw[1].nodes, fld.zw[2].nodes)
		surf = hcat(ff...)
		if showMesh
			for e = 1 : msh.N
				verts = hcat(getfield.(msh.nodes[msh.elmts[e].nodes],:xyz)...)
				plot(vcat(verts[1,:],verts[1,1]),vcat(verts[2,:],verts[2,1]), c="indigo")
			end
		end
		ax.scatter(t[1,:],t[2,:], surf, cmap="viridis", c = surf, marker = ".")
		ax.view_init(azim=0, elev=90)
		show(block=true)
	end

	function write_vtk_centres(data::DATAFIELD, dg::vecD, fi::vecD, outname::String)
		celltype = mesh_type(data.mesh) == 4 ? VTKCellTypes.VTK_QUAD : VTKCellTypes.VTK_TRIANGLE
	    cells    = Array{MeshCell}(undef,data.mesh.N)
	    nodes    = zeros(2,length(data.mesh.nodes))
		for j = 1:data.mesh.N
	        cells[j] = MeshCell(celltype, data.mesh.elmts[j].nodes)
		end
		for j = 1:length(data.mesh.nodes)
			nodes[:,j] = data.mesh.nodes[j].xyz[1:2]
		end

		vtk = vtk_grid(outname, nodes, cells)
		vtk["u_DG",VTKCellData()] = dg
		vtk["u_SIAC",VTKCellData()] = fi

		outfile = vtk_save(vtk)
	end

	function write_vtk(data::DATAFIELD, fname::String)
		celltype = mesh_type(data.mesh) == 4 ? VTKCellTypes.VTK_QUAD : VTKCellTypes.VTK_TRIANGLE
	    cells    = Array{MeshCell}(undef,data.mesh.N)
	    nodes    = zeros(2,length(data.mesh.nodes))
		for j = 1:data.N
	        cells[j] = MeshCell(celltype, data.mesh.elmts[j].nodes)
		end
		for j = 1:length(data.mesh.nodes)
			nodes[:,j] = data.mesh.nodes[j].xyz[1:2]
		end
		let outfile = vtk_grid(fname, nodes, cells, compress=3) do vtk
		    # Add some point and cell data.
			for j = 1 : size(nodes,2)
				vtk[data.field[j].name] = data.field[j].f
			end
		end
		end
	end
