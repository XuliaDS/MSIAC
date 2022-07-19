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

export element_centre, element_coords, local_coords, standard_region,standard_coords, mesh_local_to_global,
        mesh_type,create_2D_mesh,load_mesh, characteristic_size


export plot_mesh



characteristic_size(mesh::MESH) = maximum(getfield.(mesh.elmts, :area))

"""
Creates a 2D mesh file
Parameters:
===========
type: 3 = triangles, 4 = quads
nX: number of elements in the x direction
nY: number of elements in the y direction
a: domain boundary for x
b: domain boundary for y
periodic: if true, tags periodic vertices and assigns the map to its period
name: name of the mesh file

Optional:
=========
yFirst: indeces go vertically (default, we run x first)
pert: perturbate the uniform grid
structured: true or false
"""
function create_2D_mesh(type::Int, nX::Int, nY::Int, a::vecD, b::vecD, periodic::Bool, name::String ;
	                    yFirst = false , pert = [0.0,0.0], structured = true)
						if type==3 && structured == true && iszero(pert)
							return create_2D_tri_mesh(nX, nY, yFirst, a, b, periodic, name)
						end

	points = zeros((nX+1)*(nY+1),2)
	dx     = (a[2]-a[1])/nX
	dy     = (b[2]-b[1])/nY
	row    = 0
	t      = yFirst ? [b[1] .+ dy*(j-1) for j = 1:nY+1] : [a[1] .+ dx*(j-1) for j = 1:nX+1]
	n      = yFirst ? [nX+1,nY+1] : [nY+1,nX+1]

	epsX   = pert[1] == 0 ? zeros(n[2]) : rand(Uniform(-dx * pert[1],dx * pert[1]),n[2])
	epsY   = pert[2] == 0 ? zeros(n[2]) : rand(Uniform(-dy * pert[2],dy * pert[2]),n[2])
	idx = 1
	for j = 1:n[1]
		bp = [1:n[2];]
		ip = [2:n[2]-1;]
		if !structured
			(pert[1] > 0.0 ) && (epsX = rand(Uniform(-dx * pert[1],dx * pert[1]),n[2]))
			(pert[2] > 0.0 ) && (epsY = rand(Uniform(-dy * pert[2],dy * pert[2]),n[2]))
		end
		if yFirst
			points[row.+bp,1] .= a[1] .+ dx*(j-1)
			points[row.+bp,2] .= t

			points[row.+ip,2] .+= epsY[ip]
			if j != 1 && j != n[1]
				if structured
					points[row.+bp,1] .+= epsX[idx]
					idx += 1
				else
					points[row.+bp,1] += epsX
				end
			end
		else
			points[row.+bp,1] .= t
			points[row.+bp,2] .= b[1] .+ dy*(j-1)

			points[row.+ip,1] .+= epsX[ip]
			if j != 1 && j != n[1]
				if structured
					points[row.+bp,2] .+= epsY[idx]
					idx += 1
				else
					points[row.+bp,2] += epsY
				end
			end
		end
		row += n[2]
	end
	open(name, "w") do io
		println(io,"\$Nodes")
	    println(io,(nX+1) * (nY+1))
		ji = 0
		for i = 1:nY + 1
			for j = 1:nX+1
				ti = 1
				ji += 1
				if i == 1 || i == nY + 1
					ti = (j == 1 || j == nX + 1) ? -1 : 0
				elseif j == 1 || j == nX + 1
					ti = 0
				end
				println(io,ji," ", ti," ",points[ji,1]," ",points[ji,2]," ",0.0)
			end
		end
		println(io,"\$EndNodes")
		## ELEMENT MAP
		println(io,"\$Elements")
		if type == 3
			tess = basic_triangulation(points, [1:length(points[:,1]);])
			sm   = size(tess,1)
			println(io,sm)
			[ println(io,i," ", 3," ", tess[i][1]," ", tess[i][2]," ", tess[i][3]) for i = 1:sm]

			#tess = delaunay(points)
			#sm = size(tess.simplices)
			#println(io,sm[1])
			#for i = 1:sm[1]
			#	tri = tess.simplices[i,:]
			#	println(io,i," ", 3," ",tri[1]," ",tri[2]," ",tri[3])
			#end
		else
			if yFirst
		        println(io,nX * nY)
				q  = [0 nY+1 nY+2 1]
				for i = 1:nX
					[println(io,(i-1) * nY +j," ", 4," ",j+q[1]," ",j+q[2]," ",j+q[3]," ",j+q[4]," ") for j = 1:nY]
					q .+= nY+1
		        end
			else
				println(io,nX * nY)
				q  = [0 1 nX+2 nX+1]
				for i = 1:nY
					[println(io,(i-1) * nX +j," ", 4," ",j+q[1]," ",j+q[2]," ",j+q[3]," ",j+q[4]," ") for j = 1:nX]
					q .+= nX+1
		        end
			end
		end
		println(io,"\$EndElements")
		println(io,"\$Boundaries")
		println(io,4)
		print(io,"1 1 ")
		n = yFirst ? nY : nX
        [print(io, i," ") for i =1:n+1]
		print(io,"\n")

		print(io,"2 1 ")
        [print(io,i," ") for i =(nX + 1) * (nY + 1) - n:(nX + 1) * (nY + 1)]
		print(io,"\n")

		print(io,"3 1 ")
        [print(io, i," ") for i =1:n+1:(nX +1) * (nY+1)]
		print(io,"\n")

		print(io,"4 1 ")
        [print(io, i," ") for i =n+1:n+1:(nX +1) * (nY+1)]
		print(io,"\n")
		println(io,"\$EndBoundaries")
		println(io,"\$BCs")
		if periodic
			println(io,"2")
			println(io,"1 -1 1 2")
			println(io,"2 -1 3 4")
		else
			println(io,"4")
			println(io,"1 1 ")
			println(io,"2 1 ")
			println(io,"3 1 ")
			println(io,"4 1 ")
		end
		println(io,"\$EndBCs")
    end
	@info " Written in $name"
    return name
end

function create_2D_tri_mesh(nX::Int, nY::Int, orient::Bool, a::vecD, b::vecD, periodic::Bool, name::String)
    open(name, "w") do io
        println(io, "\$Nodes")
        println(io, (nX + 1) * (nY + 1))
        hx = [a[1] + j * (a[2] - a[1]) / nX for j = 0:nX]
        hy = [b[1] + j * (b[2] - b[1]) / nY for j = 0:nY]
        for i = 1:nY+1
            for j = 1:nX+1
                ti = 1
                if i == 1 || i == nY + 1
                    ti = (j == 1 || j == nX + 1) ? -1 : 0
                elseif j == 1 || j == nX + 1
                    ti = 0
                end
                println(io, (i - 1) * (nX + 1) + j, " ", ti, " ", hx[j], " ", hy[i], " ", 0.0)
            end
        end
        println(io, "\$EndNodes")
        println(io, "\$Elements")
        println(io, 2 * nX * nY)
        #Decide if diagonals have angle pi/4 or 3pi/4
        if orient
            q1 = [0 1 nX + 2]
            q2 = [0 nX + 2 nX + 1]
        else
            q2 = [1 nX + 1 0]
            q1 = [1 nX + 2 nX + 1]
        end
        for i = 1:nY
            for j = 1:nX
                println(io, 2 * (i - 1) * nX + 2 * j - 1, " ", 3, " ", j + q2[3], " ", j + q2[1], " ", j + q2[2], " ")
                println(io, 2 * (i - 1) * nX + 2 * j, " ", 3, " ", j + q1[1], " ", j + q1[2], " ", j + q1[3], " ")
            end
            q1 .+= nX + 1
            q2 .+= nX + 1
        end
        println(io, "\$EndElements")
        println(io, "\$Boundaries")

        #South
        println(io, 4)
        print(io, "1 1 ")
        [print(io, i, " ") for i = 1:nX+1]
        print(io, "\n")
        #North
        print(io, "2 1 ")
        [print(io, i, " ") for i = (nX+1)*(nY+1)-nX:(nX+1)*(nY+1)]
        print(io, "\n")
        #West
        print(io, "3 1 ")
        [print(io, i, " ") for i = 1:nX+1:(nX+1)*(nY+1)]
        print(io, "\n")
        #East
        print(io, "4 1 ")
        [print(io, i, " ") for i = nX+1:nX+1:(nX+1)*(nY+1)]
        print(io, "\n")
        println(io, "\$EndBoundaries")
        println(io, "\$BCs")
        if periodic
            println(io, "2")
            println(io, "1 -1 1 2")
            println(io, "2 -1 3 4")
        else
            println(io, "4")
            println(io, "1 1 ")
            println(io, "2 1 ")
            println(io, "3 1 ")
            println(io, "4 1 ")
        end
        println(io, "\$EndBCs")
    end
    @info "Written in $name"
    return name
end

function load_mesh(fileM::String; ordered = true::Bool, structured=false::Bool)
	 data = readdlm(fileM)
	 idN  = 1
	 N,nInfo,nType,aux = find_data("\$Nodes", data)
	 ordered && (idN += aux)
	(N == nothing) && (error("Nodes not found in mesh file\n") && return nothing)
	 E,eInfo, eType,aux = find_data("\$Elements", data[idN:end,:])
	 ordered && (idN += aux)
	(E == nothing) && (error("Elements not found in mesh file\n") && return nothing)
	L,lInfo, lType,aux = find_data("\$Boundaries",data[idN:end,:])
	ordered && (idN += aux)
	(L == nothing) && (error("Boundaries not found in mesh file\n") && return nothing)
	B,bInfo, bType,aux = find_data("\$BCs",data[idN:end,:])
	ordered && (idN += aux)
	(B == nothing) && (error("BCs not found in mesh file\n") && return nothing)

	# create a boundary map -> periodic vertex index
	pDict    = Dict{Array{Int},Array{Int}}()
	periodic = false
	bMap    = Vector{Vector{Int}}(lInfo)
	for b = 1:B[1]
		bType[b] != -1 && continue
		periodic = true
		for j = 1:length(lInfo[bInfo[b][1]])
			v1    = lInfo[bInfo[b][1]][j]
			v2    = lInfo[bInfo[b][2]][j]
			pDict = add_key(pDict, [v1], v2)
			pDict = add_key(pDict, [v2], v1)
			j == length(lInfo[bInfo[b][1]]) && continue
			b1 = sort(lInfo[bInfo[b][1]][j:j+1]) ; b2 = sort(lInfo[bInfo[b][2]][j:j+1])
			pDict = add_key(pDict, b1, b2, "val")
			pDict = add_key(pDict, b2, b1, "val")
		end
	end

	# NOW SET MESH MAP: ADJACENT ELEMENTS, NODES ETC
	vDict = Dict{Int, Vector{Int}}()
	eDict = Dict{Vector{Int},Vector{Int}}()
	for j = 1:E[1]
		piv = vcat([1:eType[j];],1)
		for f = 1:length(eInfo[j])
			vp = [eInfo[j][piv[f]], eInfo[j][piv[f+1]]]
			vDict = add_key(vDict, vp[1],j)
			eDict = add_key(eDict, sort(vp), j)
		end
	end

	m_nodes = Array{NODE, 1}(undef, N[1])
	cXY     = zeros(3) ; cX = zeros(3) ; cY = zeros(3)
	for j = 1 : N[1]
		vper = [j]
		if haskey(pDict, [j])
			aux = pDict[[j]]
			[vper = vcat(vper, pDict[[k]]) for k in aux]
			vper = vcat(vper,aux)
		end
		m_nodes[j] = NODE(nInfo[j],unique!(vDict[j]),unique!(vper))
		cj = nInfo[j]
		if j == 1
			cXY .= cj ; cY .= cj ; cX .= cj
		else
			all(cj.< cXY) && (cXY.= cj)
			cj[1] > cX[1] && (cX[1] = cj[1])
			cj[2] > cY[2] && (cY[2] = cj[2])
		end
	end
	axis = hcat([(cX .- cXY) ./ norm( cX .- cXY),
	             (cY .- cXY) ./ norm( cY .- cXY)]...)'
	m_elmts = Array{ELEMENT, 1}(undef, E[1])
	nbr = zeros(Int,(eType[1],2))
	cE  = vcat([1:eType[1];],1)
	for j = 1:E[1]
		nj   = eInfo[j]
		vj   = [m_nodes[nj[cE[k+1]]].xyz .- m_nodes[nj[cE[k]]].xyz for k = 1:eType[j]]
		area = norm(cross(vj[1], vj[2]))
		eType[j] == 3 && (area *= 0.5)
		for k = 1 : eType[j]
			nj1 = sort(nj[cE[k:k+1]])
			ele = findfirst(x->x != j, eDict[nj1])
			nbr[k,1] = ele != nothing ? eDict[nj1][ele] : haskey(pDict,nj1) ? eDict[pDict[nj1]][1] : 0
			nbr[k,2] = ele != nothing ? 1 : haskey(pDict,nj1) ? -1 : 0
		end
		m_elmts[j] = ELEMENT(eInfo[j], deepcopy(nbr), eType[j],norm.(vj), area)
	end
	return MESH(m_nodes, m_elmts, bMap, periodic, structured, axis, E[1], fileM)
end

mesh_type(mesh::MESH) = length(mesh.elmts[1].nodes)

"""
Returns the physical coordinates of the element average centre
Parameters:
===========
data: mesh & field data structure (twoD_mesh_field)
eID: element id
"""
function element_centre(mesh::MESH, eID::Int, dim = 3)
	xyz = element_coords(mesh,eID)
	return [sum(xyz[j,:]) / length(xyz[j,:]) for j = 1:dim]
end


"""
Returns the coordinates of the element vertices (ordered)
"""
element_coords(data::MESH,e::Int, dim=3) = hcat(getfield.(data.nodes[data.elmts[e].nodes],:xyz)...)[1:dim,:]


"""
Given a set of quadrature points and a quadrilateral, computes the
2D quadrature region
Params:
-------
z: quadrature points
v: element vertices
Returns:
--------
Physical coordinates and Jacobian
"""
function standard_region(z1::vecD,z2::vecD, v::matD)
	n   = size(v,2)
	i1,i2 = tensor_index(length(z1),length(z2))
	mz1   = 0.5 .*(1.0 .- z1) ; pz1 = 0.5 .*(1.0 .+ z1)
	mz2   = 0.5 .*(1.0 .- z2) ; pz2 = 0.5 .*(1.0 .+ z2)
	if n == 4
		xy = hcat([v[j,1] .* mz1[i1] .* mz2[i2] .+ v[j,2] .* pz1[i1] .* mz2[i2] .+
			       v[j,4] .* mz1[i1] .* pz2[i2] .+ v[j,3] .* pz1[i1] .* pz2[i2] for j = 1:2]...)'

		# JACOBIAN
		dz1 = [mz2[i2] .* (v[j,2].-v[j,1]) .+ pz2[i2] .* (v[j,3].-v[j,4]) for j = 1:2]
		dz2 = [mz1[i1] .* (v[j,4].-v[j,1]) .+ pz1[i1] .* (v[j,3].-v[j,2]) for j = 1:2]
		return xy, 0.25 .*(dz1[1] .* dz2[2] .- dz2[1] .* dz1[2])
	else

		xi1 = 0.5 .*(1.0 .+ z1[i1]).*(1.0 .- z2[i2]).-1.0
		xi2 = z2[i2]
		xy  = hcat(0.5 .* [v[j,1].*(.-xi2 .-xi1).+v[j,2].*(1.0.+xi1).+ v[j,3].*(1.0.+xi2) for j = 1:2]...)'

		j1 = [0.5 .* (v[j,2].-v[j,1]) for j = 1:2]
		j2 = [0.5 .* (v[j,3].-v[j,1]) for j = 1:2]
		return xy, 0.5 .* (1.0 .-z2[i2]).*(j1[1] .* j2[2] .- j2[1] .* j1[2])
		xy = hcat([v[j,1] .* mz1[i1] .* mz2[i2] .+ v[j,2] .* pz1[i1] .* mz2[i2] .+
			       v[j,3] .* pz2[i2] for j = 1:2]...)'

		dz1 = [0.5 .* mz2[i2] .* (v[j,2].-v[j,1]) for j = 1:2]
		dz2 = [0.5 .*(v[j,3] .- mz1[i1] .* v[j,1] .- pz1[i1] .* v[j,2])  for j = 1:2]
		return xy, (dz1[1] .* dz2[2] .- dz2[1] .* dz1[2])
	end
end


"""
Given a value in reference coordinates, computes the physical map
Params:
-------
z: quadrature points
v: element vertices
Returns:
--------
Physical coordinates
"""
function standard_coords(z::vecD, v::matD)
	d,n = size(v)
	mz  = 0.5 .*(1.0 .- z) ; pz = 0.5 .*(1.0 .+ z)
	if n == 4
		return v[:,1].*mz[1].*mz[2] .+ v[:,2].*pz[1].*mz[2] .+
		       v[:,4].*mz[1].*pz[2] .+ v[:,3].*pz[1].*pz[2]
	else
		return v[:,1].*mz[1].*mz[2] .+ v[:,2].*pz[1].*mz[2] .+
		       v[:,3].*pz[2]
	end
end


standard_coords_by_element(data::MESH, z::vecD, e::Int) = standard_coords(z, element_coords(data,ele))

#= this functions are used to evaluate analytic
functions: gets the corresponding parametric value
=#
function local_to_global(z1::vecD, z2::vecD, v::matD)
	d, type = size(v)
	i1,i2   = tensor_index(length(z1),length(z2))
	l2      = 0.5 .* (z2[i2] .+ 1.0)
	h       = [ norm(v[:,2] .- v[:,1]), norm(v[:,type] .- v[:,1]) ]  # 1 hori 2 vert
	t = zeros(2, length(i1))
	if type == 4
		l1 = 0.5 .* (z1[i1] .+ 1.0) ;
	else
		aux = 0.5 .*(1.0 .+ z1[i1]).*(1.0 .- z2[i2]) .- 1.0
		l1  = 0.5 .* (aux .+ 1.0) ;
	end
	t[1,:] = v[1,1].+l1.*h[2]
	t[2,:] = v[2,1].+l2.*h[1]
	return t
end

function mesh_local_to_global(data::MESH, zw1::vecD,zw2::vecD)
	q = length(zw1) * length(zw2)
	t = zeros(2,data.N * q)
	je = 1
	for e = 1:data.N
		verts = hcat(getfield.(data.nodes[data.elmts[e].nodes],:xyz)...)
		aux,dumb = standard_region(zw1,zw2,verts)
		t[1,je:je+q-1] = aux[1,:]
		t[2,je:je+q-1] = aux[2,:]
		je += q
	end
	return t
end

function local_coords(p,v)
	d, n = size(v)
	d == 2 && (v = vcat(v, zeros(n)'))
	dv1 = v[:,2].- v[:,1] ;	dv2 = v[:,n].- v[:,1]

	np   = cross(dv1,dv2)
	p3   = length(p) == 2 ? vcat(p,0.0) : p
	o1   = cross(np,dv1); o2 = cross(np,dv2)
	z    = [2.0 .* dot(p3 .- v[:,1],o2) / dot(dv1,o2).-1.0, 2.0 .* dot(p3 .- v[:,1],o1) / dot(dv2,o1).-1.0]
	if (n == 4) return z
	elseif abs(1.0 -z[2]) < fTOL return [-1.0,z[2]]
	else return [(1 .+ z[1]) ./ (1.0 .- z[2]) .* 2.0 .- 1.0, z[2]]
	end
end



## DEBUG FUNCTION


function plot_mesh_lineSupport(data::MESH, point::vecD,foot::Vector{Matrix{Float64}}; saveFig = nothing::Union{Nothing,String}, name="")
		fig = figure(5,figsize=(10,10))
		ax  = fig.gca()
		if name != ""
			ax.set_title(name)
		end
		nc   = 7
		cmap = matplotlib.cm.inferno([0:1.0 / (nc-1):1;])
		c1 = cmap[6,:] ;c2 = cmap[1,:]; c3 = cmap[2,:]; c4 = cmap[4,:]
		for j = 1:data.N
			verts = data.elmts[j].nodes
			xyz   = hcat(getfield.(data.nodes[verts],:xyz)...)
			lv    = length(xyz[1,:])
			aux   = ele_ca(lv)
			scatter(xyz[1,:],xyz[2,:], color = "w", s=50)
			ax.plot(xyz[1,aux], xyz[2,aux], color = c1, lw = 2)
			ave = [1 / lv * sum(xyz[k,:]) for k = 1:3]
			ax.text(ave[1], ave[2],j, fontsize=12, zorder = 5, color = "black")
		end
		for j = 1:length(data.nodes)
			ax.text(data.nodes[j].xyz[1,:],data.nodes[j].xyz[2,:],j, fontsize=12, zorder = 5, color = "red")
		end
		tr = size(foot)
		for e = 1:tr[1]
			ir = foot[e]
			ca = ele_ca(length(ir[1,:]))
			ax.plot(ir[1,ca],ir[2,ca], color = c3, ls = "--", lw = 2)
			ax.fill(ir[1,ca],ir[2,ca], edgecolor="w", facecolor = c3, zorder = 0, alpha = 0.15)
	     ax.scatter(ir[1,:],ir[2,:], color = c3)
		end
		ax.scatter(point[1],point[2], color = c4, marker = "o", s = 180,zorder = 4)
		ax.scatter(point[1],point[2], color = "white", marker = "o", s = 160,zorder = 4)
		ax.scatter(point[1],point[2], color = c4, marker = "*", s = 110,zorder = 5)
		margins = ["right","top","bottom", "left"]
		[ax.spines[margins[i]].set_visible(false) for i in 1:length(margins)]
		ax.set_yticks([])
	 	ax.set_xticks([])
		if saveFig != nothing
			savefig(saveFig*".pdf", bbox_incMSIAC_hes="tight")
			@info " Written on $saveFig .pdf "
			close(fig)
		else
	    	return
			show(block=true)
		end
	end
	function plot_mesh(data::Union{DATAFIELD,MESH}; blockPlot = true::Bool, noIds = false::Bool, savePlot = " "::String)
		if typeof(data) == DATAFIELD
			data = data.mesh
		end
		fig = figure(1)
		ax   = fig.gca()
		nc   = 7
		cmap = matplotlib.cm.inferno([0:1.0 / (nc-1):1;])
		for j = 1:data.N
			verts = data.elmts[j].nodes
			xyz   = hcat(getfield.(data.nodes[verts],:xyz)...)
			lv    = length(xyz[1,:])
			aux   = ele_ca(lv)
			ax.scatter(xyz[1,:],xyz[2,:], color = "w", s=50)
			ax.plot(xyz[1,aux], xyz[2,aux], color = cmap[3,:], lw = 0.75)
			ave = [1 / lv * sum(xyz[k,:]) for k = 1:3]
			noIds  && continue
			ax.text(ave[1], ave[2],string(j), fontsize=12)
			[ax.text(xyz[1,k], xyz[2,k], string(verts[k])) for k = 1:lv]
		end
		margins = ["right","top","bottom", "left"]
		[ax.spines[margins[i]].set_visible(false) for i in 1:length(margins)]
		plt.yticks([])
	    plt.xticks([])
		if savePlot != " "
			savefig(savePlot*".pdf", bbox_inches="tight")
			@info " Written on $savePlot .pdf "
			close(fig)
		else
			show(block=blockPlot)
		end
	end
