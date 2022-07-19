#=
*
*      MSIAC: a Julia module designed to post-process data with SIAC filters
*                                Version 1.0
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
=#
__precompile__()
#using Revise
module MSIAC

using LinearAlgebra, Einsum,DelimitedFiles, Triangle, Distributions, Printf,CircularArrays, Jacobi, SharedArrays, Distributed
import BSplines, FastGaussQuadrature, Jacobi, Polynomials
const JP   = Jacobi
const GQ   = FastGaussQuadrature


const fTOL = 1.e-12
const dbl  = Float64
const vecD = Vector{Float64}
const vecI = Vector{Int}
const matD = Matrix{Float64}
const matI = Matrix{Int}

ele_ca(type::Int) = repeat([1:type;],2)

export  get_eL2, test_filter, GZW, in_polygon, evaluate_point, polynomial_degree, MESH, FIELD, DATAFIELD,
        tensor_index, filter_convolution, filter_data

function tensor_index(n1::Int,n2::Int, reverse = false::Bool)
	x = vcat([fill(j,n2) for j = 1:n1]...)
	y = repeat([1:n2;], n1)
	if reverse
		return y,x
	else
		return x,y
	end
end

"""
Structure for quadrature
------------------------
nodes :
	quadrature points
weights:

eQT:
	type: Legendre, Lobatto, etc
"""
struct GZW
   type    :: String
   Q       :: Int
   nodes   :: vecD
   weights :: vecD
   function GZW(tQ::String, Q::Int)
	   if tQ == "Legendre" || tQ == "legendre" || tQ == "LEGENDRE"
	       z,w = GQ.gausslegendre(Q)
	   elseif tQ == "Lobatto" || tQ == "lobatto" || tQ == "LOBATTO"
	       z,w = GQ.gausslobatto(Q)
	   elseif tQ == "Radau" || tQ == "radau" || tQ == "RADAU"
		   z,w = GQ.gaussradau(Q)
	   else
	       print(" Quadrature type ", tQ," not supported. Add it ! \n")
	       return nothing
	   end
	   return new(tQ, Q, z, w)
   end
end

"""
holds footprint information: element reference coordinates (z)
kernel coordinates (t) and sorts by elements (e)
"""
mutable struct BP
	z :: vecD
	t :: vecD
	e :: vecI
	v :: Int
	BP(z::vecD,t::vecD,e::vecI,v :: Int) = new(deepcopy(z), deepcopy(t),deepcopy(e), v)
end

struct POINT
    x::Float64
    y::Float64
	POINT(x::vecD) = new(x[1],x[2])
	POINT(x1::dbl,x2::dbl) = new(x1,x2)
end
const SEG = Union{Tuple{POINT,POINT},Tuple{vecD,vecD}}
Base.:+(p::POINT,q::POINT) = (p.x+q.x,p.y+q.y)
Base.:-(p::POINT,q::POINT) = (p.x-q.x,p.y-q.y)
Base.:+(p::POINT,q::vecD)  = (p.x+q[1],p.y+q[2])
Base.:-(p::POINT,q::vecD)  = (p.x+q[1],p.y+q[1])

struct EXPANSION
	type     :: String
	degree   :: Int
	J        :: Vector{Polynomials.Polynomial}
	mass     :: matD
	collapse :: Bool
	m1       :: vecI
	m2       :: vecI
	function EXPANSION(type::String,p::Int, collapse::Bool)
		if type == "Pk" || type == "pk" || type == "PK"
			m1 = vcat([fill(i,p+1-i) for i = 0:p]...)
			m2 = vcat([[0:p-i;] for i = 0:p]...)
			return new(type,p,Vector{Polynomials.Polynomial}(undef,1),ones(1,1),collapse,m1,m2)
		else
			m1,m2 = tensor_index(p+1,p+1)
			if type == "legendre" || type == "LEGENDRE" || type == "Legendre"
				J  = [JP.poly_legendre(i) for i = 0:p]
			else
				PJ = [JP.poly_jacobi(i,1.0,1.0) for i = 0:p-2]
				J   = vcat(Polynomials.Polynomial([0.5,.-0.5]),
				          [Polynomials.Polynomial([0.25,0.0,.-0.25]) * PJ[i] for i = 1:p-1],
					       Polynomials.Polynomial([0.5,0.5]))
		    end
			# Get mass matrix
			q1    = GZW("legendre", p + 1)
			q2    = GZW("legendre", p + 1)
			t1,t2 = tensor_index(q1.Q , q2.Q)
			pl2   = zeros((p+1)^2,(p+1)^2)
			n1    = q1.nodes[t1] ;	n2    = q2.nodes[t2]
			we    = q1.weights[t1] .* q2.weights[t2]
			for i = 1:length(m1)
				p1 = J[m1[i]].(n1) .* J[m2[i]].(n2)
				for j = 1:length(m1)
					p2 = J[m1[j]].(n1) .* J[m2[j]].(n2)
					pl2[i,j] = sum(p1 .* p2 .* we)
				end
			end
			mij  = inv(pl2)
			return new(type,p,J,mij,collapse,m1,m2)
		end
	end
end

"""
Modal reconstruction based on expansion (basis) and f sampled at zw12 points
Parameters:
===========
basis: a EXPANSION structure (basis, mass matrix,etc)
f: the projection function at a particular element size length(zw1) * length(zw2)
zw1: points sampled along 1st dir
zw2: points sampled along 2nd dir

Returns:
========
modes matrix with dimension [data.degree + 1, #elements]
"""
function modal_values(basis::EXPANSION, f::matD, zw::Vector{GZW}, reverse=false::Bool)
	# set mass matrix
	Q,n = size(f)
	i1,i2 = tensor_index(zw[1].Q , zw[2].Q, reverse)
	if (Q != length(i1))
		@error " Inside modes_f :: MISMATCH BETWEEN WEIGHTS AND FIELD $Q != $(zw[1].Q * zw[2].Q)"
		return nothing
	end
	w12  = zw[1].weights[i1] .* zw[2].weights[i2]
	phi  = eval_basis(basis,zw[1].nodes[i1],zw[2].nodes[i2])
	if basis.type == "Pk" || basis.type == "pk" || basis.type == "Pk"
		return [sum(f[:, e] .* w12 .* phi[:, j] .* 0.5 .* (1 .- zw[2].nodes[i2])) for j = 1:length(basis.m1), e = 1:n]
	else
		return hcat([basis.mass * [sum(f[:,e] .* phi[:,j] .* w12) for j = 1 : length(phi[1,:])] for e = 1:n]...)
	end
end

function phys_values(basis::EXPANSION, modes::matD, zw::Vector{GZW}, reverse=false::Bool)
	i1, i2 = tensor_index(zw[1].Q, zw[2].Q,reverse)
	P      = eval_basis(basis,zw[1].nodes[i1],zw[2].nodes[i2])
	proj   = zeros(length(i1),size(modes,2))
	@einsum proj[k,e] = modes[j,e] * P[k,j]
	return proj
end


mutable struct FIELD
	zwC     :: Vector{GZW}
	zw      :: Vector{GZW}
	reverse :: Bool
	f       :: matD
	modes   :: matD
	name    :: String
	basis   :: EXPANSION
	function FIELD(type::String, p::Int, collapse::Bool; reverse=false::Bool,zw = GZW[]::Vector{GZW},
		            f = zeros(0,0)::matD, m = zeros(0,0)::matD,n = ""::String)
				   basis = EXPANSION(type,p,collapse)
				   qt = ["legendre","legendre"] ; Q = p + 2
		   		   collapse && (qt[2] = "Radau")
		   		   isempty(zw) && (zw = [GZW(qt[1],p+1),GZW(qt[2],p+1)])
				   ff = !isempty(f) ? f : !isempty(m) ? phys_values(basis,m,zw,reverse)  : zeros(0,0)
				   mm = !isempty(m) ? m : !isempty(f) ? modal_values(basis,f,zw,reverse) : zeros(0,0)
				   return new(deepcopy(zw),deepcopy(zw),reverse,ff,mm,n,basis)
	end
end

struct NODE
	xyz  :: vecD
	ele  :: vecI
	type :: vecI
end

struct ELEMENT
	nodes :: vecI
	neigh :: matI
	type  :: Int
	h     :: vecD
	area  :: dbl
end

struct MESH
	nodes      :: Vector{NODE}
	elmts      :: Vector{ELEMENT}
	bMap       :: Vector{Vector{Int}}
	periodic   :: Bool
	structured :: Bool
	frame      :: matD
	N          :: Int
	fileM      :: String
end

mutable struct DATAFIELD
	fields :: Vector{FIELD}
	fdict  :: Dict{String,Int}
	mesh   :: MESH
	function DATAFIELD(f::Vector{FIELD}, m::MESH)
		fdict = Dict("" => 1)
		[fdict[f[k].name] = k for k = 1: length(f)]
		return new(f,fdict,m)
	end
end

polynomial_degree(data::DATAFIELD) = data.fields[data.fdict[""]].basis.degree

function find_data(name::String, data::Matrix{Any})
	i1 = findfirst(x->x == name, data)[1]
	i2 = findfirst(x->x == "\$End"*strip(name,['\$']), data)[1]
	if i1 == nothing || i2 == nothing
		@error " !!!!! MISSING KEYWORD $name IN DATA MATRIX \n"
		return nothing
	end
	info     = data[i1 + 1,:]
	dType    = zeros(Int,info[1])
	datatype = typeof(data[i1+2,3])
	dVals    = Vector{datatype}[[] for m = 1:info[1]]
	for j = i1 + 2:i2-1
		j_st = findfirst(x->x == "",data[j,:])
		dl   =  j_st == nothing ? length(data[j,:]) : j_st - 1 #- length(x)
		dVals[data[j,1]] = datatype.(data[j,3:dl])  # data not necessarily in order
		dType[data[j,1]] = data[j,2]     # data not necessarily in order
	end
	info, dVals, dType, i2
end

function add_key(d,k, val, type = "vector")
	if haskey(d,k)
		push!(d[k],val)
	elseif type == "vector"
		d[k] = [val]
	else
		d[k] = val
	end
	return d
end

include("kernel.jl")
include("mesh.jl")
include("fields.jl")
include("structured_triangles_scripts.jl")

## Geometry functions

cross2D(v1::vecD, v2::vecD) = v1[1] * v2[2] - v1[2] * v2[1]
ndot(v1::vecD,v2::vecD) = dot(v1./norm(v1),v2 ./ norm(v2))

in_segment(p::vecD, s::SEG, zeroTOL = fTOL) = abs(norm(p-s[1]) + norm(s[2]-p) - norm(s[2]-s[1])) < zeroTOL
function edges_disjoint(ls1::SEG, ls2::SEG; witness::Bool=false, zeroTOL::Real=fTOL)

		r = ls1[2] - ls1[1]
		s = ls2[2] - ls2[1]
	    # first line segment is a point
		if norm(r) < zeroTOL
	        empty_intersection = !in_segment(ls1[2],ls2)
	        if witness
	            return (empty_intersection, empty_intersection ? dbl[] : ls1[2])
	        else
	            return empty_intersection
	        end
	    end
		if norm(s) < zeroTOL
	        # second line segment is a point
	        empty_intersection = !in_segment(ls2[2],ls1)
	        if witness
	            return (empty_intersection, empty_intersection ? dbl[] : ls2[2])
	        else
	            return empty_intersection
	        end
	    end

	    p1p2 = ls2[1] - ls1[1]
	    u_numerator   = cross2D(p1p2, r)
	    u_denominator = cross2D(r, s)

	    if abs(u_denominator) < zeroTOL
			empty_intersection = true
	        # line segments are parallel
	        if abs(u_numerator) < zeroTOL
	            # line segments are collinear
				if in_segment(ls1[2],ls2)
	                empty_intersection = false
	                witness && (v = ls1[2]) #: [ls1[1], ls1[2]]

				elseif in_segment(ls2[1],ls1)
	                empty_intersection = false
					witness && (v = ls2[1]) #: [ls1[1], ls1[2]]

				elseif in_segment(ls2[2],ls1)
	                empty_intersection = false
					witness && (v = ls2[2]) #: [ls1[1], ls1[2]]

	            elseif in_segment(ls1[1],ls2)
	                empty_intersection = false
					witness && (v = ls1[1]) #: [ls1[1], ls1[2]]
	            end
			end
	    else
	        # line segments are not parallel
	        u = u_numerator / u_denominator
	        if u < -zeroTOL || u > 1.0 + zeroTOL
	            empty_intersection = true
	        else
	            t = cross2D(p1p2, s) / u_denominator
	            empty_intersection = t < -zeroTOL || t > 1.0 + zeroTOL
	            if witness
	                v = ls1[1] + t * r
	            end
	        end
	    end
	    if witness
	        return (empty_intersection, empty_intersection ? dbl[] : v)
	    else
	        return empty_intersection
	    end
	end

# ** Graham scan algorithm : Sorting points ccw  **

ccw(a::POINT, b::POINT, c::POINT) = ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x))
ccwV(a::vecD, b::vecD, c::vecD) = ((b[1] - a[1])*(c[2] - a[2]) - (b[2] - a[2])*(c[1] - a[1]))

function graham_scan!(points::Vector{POINT})
	N = length(points)
	# Place the lowest point at the start of the array
    p = sortperm(points, by = item -> item.y)
	# Sort all other points according to angle with that point
    q = sortperm(points[p[2:end]], by = item -> atan(item.y - points[p[1]].y,item.x - points[p[1]].x))

	p[2:end] = p[q.+1]
	#@info " CONVEX POLY INITAL PIV $p "
	# M will be the point on the hull
    M = 2
	for i = 1:N
		do_while = false
        while (ccw(points[p[M-1]], points[p[M]], points[p[i]]) < 0.0)
			do_while = true
            if (M > 2)
				#@info " REDUCE BY 1 $M "
                M -= 1
            # All points are collinear
            elseif (i == N)
				#@info " ALL BPS COLLINEAR "
                break
            else
				#@info " NEXT "
                i += 1
            end
        end
        # ccw point found, updating hull and swapping points
        M += 1
		if do_while
			#@info " M NOW $M i $i swap --- PIV $p"
			p[i],p[M] = p[M], p[i]
		end
		M == N && break
    end
	return p
end

convex_polygon(v::Vector{Vector{dbl}}) = graham_scan!([POINT(v[j]) for j = 1:size(v,1)])
# *************************************************
function rayintersectseg(p::POINT, edge::SEG)
	a, b = edge
	if a.y > b.y
		a, b = b, a
	end
	#@info " POINT $p"
	px = p.x ∈ (a.x,b.x) ? p.x + eps(p.x) : p.x# == a.x ? p.x + eps(p.x) : p.x == b.x ? p.x - eps(p.x) : p.x
	py = p.y ∈ (a.y,b.y) ? p.y + eps(p.y) : p.y # == a.y ? p.y + eps(p.y) : p.y == b.y ? p.y - eps(p.y) : p.y
	p = POINT([px,py])
	if (p.y > b.y || p.y < a.y) || (p.x > max(a.x, b.x))
		return false
	end

	if p.x < min(a.x, b.x)
		true
	else
		mred = (b.y - a.y) / (b.x - a.x)
		mblu = (p.y - a.y) / (p.x - a.x)
		return mblu ≥ mred
	end
end

function in_polygon(p::vecD, verts::matD ; zeroTOL = fTOL::dbl)
	nv   = size(verts,2)
	ca   = CircularArray([1:nv;])
	hull = [0,0]
	for j = 1:nv
		norm(p-verts[1:2,j]) < zeroTOL && return [1,j]
		in_segment(p, (verts[1:2,ca[j]], verts[1:2,ca[j+1]]), zeroTOL) && (hull = [0,j])
	end
	!iszero(hull) && return hull
	xl = [minimum(verts[1,:]),maximum(verts[1,:])]
	yl = [minimum(verts[2,:]),maximum(verts[2,:])]
	(p[1] < xl[1] || p[1] > xl[2] || p[2] < yl[1] || p[2] > yl[2]) && return -1
	pt = POINT(p)
	rays = [rayintersectseg(pt, (POINT(verts[:,j]),POINT(verts[:,ca[j+1]])) ) for j = 1 : nv]
	return isodd(count(rays)) ? 0 : -1
end
##############

function periodic_point!(data::MESH,v::vecI, pC::BP, e2, q::vecD)
	d   = 2
	lpq = (pC.z,q)
	xy  = (data.nodes[v[1]].xyz[1:d],data.nodes[v[2]].xyz[1:d])
	res = edges_disjoint(lpq, xy ; witness = true)
	res[1] == true && @error " periodic point :: segment pq $lpq doesn't cross vertices $v"
	t   = norm(res[2] .- xy[1]) / norm(xy[2] .- xy[1])  #WARNING THIS IS ASSUMING STRAIGHT ELEMENTS
	v1  = findfirst(x->x ∈ data.nodes[v[1]].type, data.elmts[e2].nodes)
	v2  = findfirst(x->x ∈ data.nodes[v[2]].type, data.elmts[e2].nodes)
	if v1 == nothing || v2 == nothing
		@error " This is not a periodic edge !!!! "
		return
	end
	pC.z = data.nodes[data.elmts[e2].nodes[v1]].xyz[1:d] .+
	       t .* (data.nodes[data.elmts[e2].nodes[v2]].xyz[1:d] .- data.nodes[data.elmts[e2].nodes[v1]].xyz[1:d])
    return pC.z .+ (q[1:d].-res[2])
end

function add_break!(data::MESH, breaks::Vector{BP}, pt::BP)

	isempty(breaks) && (breaks = Vector{BP}([BP(pt.z,pt.t,pt.e,pt.v)]))
	piv = findfirst(x -> norm(x .- pt.z) < fTOL, getfield.(breaks,:z) )
	if piv == nothing
		push!(breaks, BP(pt.z,pt.t,pt.e,pt.v))
	else
		#@info " ADD element $(pt.e) to existing break at $piv "
		breaks[piv].e = unique!(vcat(breaks[piv].e,pt.e))
	end
	return breaks
end

function wall_distance!(data::MESH, breaks::Vector{BP}, pC::BP, cN::vecD,tDir::vecD)
	verts = element_coords(data,pC.e[1],2)
	cPN   = (pC.z,cN)
	qc    = ele_ca(data.elmts[pC.e[1]].type)
	for j = 1:length(verts[:,1])
		(data.elmts[pC.e[1]].neigh[j,1] > 0) && continue
		edge = (verts[:,qc[j]],verts[:,qc[j+1]])
		if abs(ndot(edge[2] .- edge[1], tDir) ) > 0.9999
			continue
		end
		res  = edges_disjoint(cPN,edge; witness = true)
		res[1] && continue
		pC.z = res[2];
		pC.v = norm(pC.z .- verts[:,qc[j]])   < fTOL ? data.elmts[pC.e[1]].nodes[qc[j]] :
		       norm(pC.z .- verts[:,qc[j+1]]) < fTOL ? data.elmts[pC.e[1]].nodes[qc[j+1]] : 0
		tACC  = norm(cPN[1] .- pC.z) ./ norm(tDir)
		pC.t += norm(cPN[1] .- pC.z) .* tDir
		return tACC, add_break!(data, breaks,pC)
	end
	return 0.0, breaks
end

function pull_vertex!(data::MESH, pC::BP, dir::vecD)
	vper  = data.nodes[pC.v].type
	veles = vcat(getfield.(data.nodes[vper],:ele)...)
	vmain = []
	[ (pC.v in data.elmts[e].nodes) && push!(vmain,e) for e in veles]

	aux = findall(x->x ==  1, data.elmts[pC.e[1]].neigh[:,2])
	o1  = findall(x->x in data.elmts[pC.e[1]].neigh[aux,1], vmain)
	aux = findall(x->x == -1, data.elmts[pC.e[1]].neigh[:,2])
	o2  = findall(x->x in data.elmts[pC.e[1]].neigh[aux,1], veles)
	eles = unique(vcat(pC.e[1],vmain[o1],vmain, veles[o2],veles))
	for j in eles
		le   = data.elmts[j]
		vs   = element_coords(data,j,2)
		aux  = findfirst(x->x ∈ vper, le.nodes)
		aux == nothing && continue
		if le.nodes[aux] != pC.v && abs(ndot(vs[:,aux] .- pC.z, dir) ) < 0.99
			continue
		end
		cPN = (vs[:,aux],vs[:,aux].+ dir)
		ca  = ele_ca(le.type)
		for k = 1 : le.type
			edge = (vs[:,ca[k]], vs[:,ca[k+1]])
			foundV = false
			if any(ca[k:k+1] .== aux)
				s1 = in_segment(edge[1],cPN)
				s2 = in_segment(edge[2],cPN)
				s1 && s2 && (foundV = true )
			else
				!edges_disjoint(cPN,edge) && (foundV = true)
			end
			if foundV
				pC.e[1] = j ; pC.z = vs[:,aux] ; pC.v = le.nodes[aux]
				return
			end
		end
	end
	pC.e[1] = 0
end

function support_fits!(data::MESH, pC::BP, dir::vecD, T::vecD;
	                   breaks = BP[]::Vector{BP})

	(length(T) == 0) && (return BP[], 0.0)
	## for knots: we will have unscaled breaks
	ndir      = dir ./ norm(dir)
	tDir      = sign(T[end]) .* ndir
	largeStep = max(10.0, abs(T[end])) * maximum(data.elmts[pC.e[1]].h) * tDir
	cN   = pC.z .+ largeStep
	tPN  = 0.0;
	cPN  = (pC.z,cN)
	eAux = -1
	# Check if pp point is in element hull: need to choose the appropriate starting element
	ca  = ele_ca(data.elmts[pC.e[1]].type)
	vs  = element_coords(data,pC.e[1],2)
	pcb = pC.e[1]
	for j = 1:length(vs[1,:])
		if norm(pC.z .- vs[:,j]) < fTOL
			pC.v = data.elmts[pC.e[1]].nodes[j]
			eAux = 1
			pull_vertex!(data, pC, largeStep)
			if pC.e[1] == 0
				pC.e[1] = pcb
				dist, breaks = wall_distance!(data, breaks,pC,cN,tDir)
				return breaks, abs(dist - T[end]) < fTOL ? 0.0 : dist-T[end]
			end
			break
		end
	end
	if eAux == -1
		ip = [0,0]
		for j = 1:length(vs[1,:])
			edge = (vs[:,ca[j]], vs[:,ca[j+1]])
			if in_segment(pC.z, edge)
				ip[1] = j
				vn1   = (edge[2].-edge[1]) / norm(edge[2] .- edge[1])
				vn2   = (cN .- pC.z) / norm(cN .- pC.z)
				if abs(dot(vn1,vn2)) > 0.99999
					ip[1] = 0
					break
				end
			elseif !edges_disjoint(cPN,edge)
				ip[2]  = j
			end
			all(ip .!= 0) && break
		end
		if ip[1] != 0 && ip[2] == 0
			if data.elmts[pC.e[1]].neigh[ip[1],2] == -1
				le = data.elmts[pC.e[1]]
				ca = ele_ca(le.type)
				per = le.nodes[ca[ip[1]:ip[1]+1]]
				cN = periodic_point!(data,per, pC, le.neigh[ip[1],1], cN)
			end
			pC.e[1] = data.elmts[pC.e[1]].neigh[ip[1],1]
		end
		(ip[1] > 0 ) && (eAux = 0)
		if pC.e[1] == 0
			pC.e[1] = pcb
			dist, breaks = wall_distance!(data, breaks, pC, cN, tDir)
			return breaks, abs(dist - T[end]) < fTOL ? 0.0 : dist-T[end]
		end
	end

	back = -1
	eP   = zeros(2)
	(eAux >= 0) && (breaks = add_break!(data, breaks,pC))
	#WARNING THIS WILL HAVE TO CHANGE WE HAVE MESHES WITH MIXED ELEMENTS
	ca = ele_ca(data.elmts[pC.e[1]].type)
	##
	for knot in T
		cN  = pC.z .+ (knot - tPN) .* dir
		(eAux >=0 && abs(knot) < fTOL) && continue # skip, break was counted at start
		do_while = true
		#find all intersections between two spline knots
		while (do_while)
			pC.v = 0; eAdj = 0; alongEdge = false
			cPN = (pC.z,cN)
			# loop around element: find intersection point
			le = data.elmts[pC.e[1]]
			lv = element_coords(data,pC.e[1],2)
			inPol  = in_polygon(cN, lv, zeroTOL = fTOL)
			inEdge = [in_segment(pC.z,(lv[:,ca[j]], lv[:,ca[j+1]])) for j = 1: le.type]
			if typeof(inPol) == vecI
				eP = cN[1:2]
				eAdj      = inPol[2] ;
				alongEdge = inEdge[inPol[2]]
				inPol[1] == 1 && (pC.v = le.nodes[inPol[2]])
			elseif inPol == -1
				for j = 1: le.type
					(back != -1 && le.neigh[j,1] == back) && continue
					edge = (lv[:,ca[j]], lv[:,ca[j+1]])
					disj,eP = edges_disjoint(cPN,edge; witness = true)
					if !disj
						(norm(eP .- pC.z) < fTOL) && continue
						pC.v = norm(edge[1]- eP) < fTOL ? le.nodes[ca[j  ]] :
							   norm(edge[2]- eP) < fTOL ? le.nodes[ca[j+1]] : 0
						eAdj = j
						break
					end
				end
			end
			iP = eAdj == 0 ? deepcopy(cN) : deepcopy(eP)
			# check that point is actually in element (can be in adjacent for unstructured meshes )
			if !alongEdge && eAdj == 0 && inPol == -1
				id = findall(inEdge .== true)
				pC.e[1] = data.elmts[pC.e[1]].neigh[id[1],1]
				breaks = add_break!(data, breaks,pC) #add break from current element
				continue
			end
			nt    = norm(pC.z .- iP)
			tACC  = nt * sign(T[end]) / norm(dir) # kernel coordinates
			pC.t += nt .* tDir
			tPN  += tACC
			pC.z .= iP
			breaks   = add_break!(data, breaks,pC) #add break from current element
			abs(tPN - T[end]) < fTOL && break
			do_while = norm(pC.z .- cN) > fTOL
			(eAdj == 0) && continue
			next = (alongEdge && pC.v <1) ? pC.e[1] : le.neigh[eAdj,1]
			back = pC.e[1]

			if pC.v > 0 # intersection point is vertex. could have wrong neighbor
				cNP = norm(pC.z .- cN)
				pull_vertex!(data, pC, largeStep)
				next = pC.e[1]
				cN .= pC.z .+ cNP.* tDir  # if we have periodic bounds, need to update coordinates from per. vertex
			elseif !alongEdge && next != 0 && le.neigh[eAdj,2] == -1
				cN = periodic_point!(data,le.nodes[ca[eAdj:eAdj+1]], pC, next, cN)
			end
			if next == 0  # we have hit a wall. Return shift for kernel
				pC.e[1] = back
				dist, breaks = wall_distance!(data, breaks, pC,cN,tDir)
				tPN += dist
				return breaks, abs(tPN - T[end]) < fTOL ? 0.0 : tPN-T[end]
			end
			pC.e[1] = next
			next != back && (breaks = add_break!(data, breaks,pC)) #add break from next element
		end
	end
	return breaks, abs(tPN - T[end]) < 1.e-14 ? 0.0 : tPN-T[end]
end




"""
Given a mesh & field data structure, filters the entire field
Parameters:
===========
data: mesh & field data structrue (FIELD)
kType: kernel type (line / tensor)

Optional:
=========
rx: number of splines in the 1st direction / line direction

lx: order of splines in the 1st direction / line direction

ry: number of splines in the 2nd direction (only useful for tensor filter)

ly: order of splines in the 2nd direction (only useful for tensor filter)

dumpFile: if true, writes an output file with the filtered data
"""
function filter_data(msh::MESH, fld::FIELD, kType::String; rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
	                 dumpFile="nofile"::String, theta="auto"::Union{String,dbl,Vector{dbl}},
					 scaling="auto"::Union{String,dbl,Vector{dbl}},inexact=false::Bool, parallel=true::Bool)

	(rx == -1) && (rx = 2 * fld.basis.degree + 1)
	(ry == -1) && (ry = 2 * fld.basis.degree + 1)
	(lx == -1) && (lx =     fld.basis.degree + 1)
	(ly == -1) && (ly =     fld.basis.degree + 1)
	kx = [rx,lx,rx,lx,1] ; ky =[ry,ly,ry,ly,1]
	Q          = max(2,ceil(Int, 0.5 * (3 * fld.basis.degree + max(lx, ly) + 1)))
	fld.zwC[1] = GZW("legendre",Q)
	fld.zwC[2] = fld.basis.collapse ? GZW("Radau",Q) : GZW("legendre",Q)

	i1,i2 = tensor_index(fld.zw[1].Q,fld.zw[2].Q,fld.reverse)
	Qt    = length(i1)
	filt  = zeros(Qt,msh.N)
	if parallel 	
		i1C = CircularArray(fld.zw[1].nodes[i1])
		i2C = CircularArray(fld.zw[2].nodes[i2])
		eC  = vcat([fill(e,Qt) for e in 1:msh.N]...)

		aux = zeros(Qt*msh.N)
		Threads.@threads for e in 1:msh.N * Qt
			aux[e] = filter_convolution(msh,fld,eC[e],[i1C[e],i2C[e]],kType,kx = kx, ky = ky,theta = theta, scaling = scaling, inexact=inexact)
		end
		filt = reshape(aux, Qt, msh.N)
			
	else
		for e in 1:msh.N
			for i in 1:Qt
				filt[i,e] = filter_convolution(msh,fld,e,[fld.zw[1].nodes[i1[i]],fld.zw[2].nodes[i2[i]]],kType,kx =kx, ky = ky,theta = theta, scaling = scaling, inexact=inexact)
			end
		end
	end

	dumpFile == "nofile" && (return filt)

	ext = split(dumpFile, ".")[end]
	ext == dumpFile && (dumpFile *=".txt")
	open(dumpFile, "w") do io
		println(io,"\$Fields")
		println(io,msh.N," ", fld.zw[1].Q," ", fld.zw[1].type," ", fld.zw[2].Q," ",fld.zw[2].type)
		for j = 1:msh.N
			print(io,j," ",1," ")
			for k = 1:fld.zw[1].Q * fld.zw[2].Q
				print(io,filt[k,j]," ")
			end
			print(io,"\n")
		end
		println(io,"\$EndFields")
	end
	print(" \n -------- WRITTEN ON FILE ", dumpFile,"\n")
	#showWarnings && @info " Finished postprocessing field with $kType kernel. Total points: $(data.N * data.zw[1].Q * data.zw[2].Q) "
	return filt
end

filter_data(data::DATAFIELD, kType::String; f=""::String,rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
	                 dumpFile="nofile"::String, theta="auto"::Union{String,dbl,Vector{dbl}},
					 scaling="auto"::Union{String,dbl,Vector{dbl}}, parallel=true::Bool, inexact=false::Bool) =
					 filter_data(data.mesh,data.fields[data.fdict[f]],kType, rx = rx, lx = lx, ry = ry, ly=ly,
					 dumpFile = dumpFile, theta=theta, scaling=scaling, parallel=parallel, inexact=inexact)

function get_element_list!(ker::Vector{KERNEL}, data::MESH, pC::BP, kdir::vecD;
	                       kT = [eSYM, ePOS]::vecI,secondAxis = false::Bool, skip_knots = false)

	eqK = length(kT) == 2 && ker[1].rC == ker[2].rC && ker[1].l == ker[2].l
	for t in kT
		update_knot_matrix!(ker[t])
		kL,kR = skip_knots ? ([minimum(ker[t].T)],[maximum(ker[t].T)]) : split_knots(ker[t].T)
		b0 = secondAxis ? add_break!(data,BP[],pC) : BP[]
		bl1, dL = support_fits!(data,deepcopy(pC), kdir, kL,breaks = b0)

		(t == eSYM && dL > 0.0 && !eqK) && continue
		bline, dR = support_fits!(data,deepcopy(pC), kdir, kR,
		                          breaks= isempty(bl1) ? b0 : reverse!(bl1))
		(t == eSYM && dR < 0.0 && !eqK) && continue
		shift = dL > 0.0 ? dL : dR < 0.0 ? dR : 0.0
		if shift != 0.0
			length(kT) == 2 && (t = ePOS)
			update_knot_matrix!(ker[t]; lX = shift)
			kL,kR = skip_knots ? ([minimum(ker[t].T)],[maximum(ker[t].T)]) : split_knots(ker[t].T)
			b0 = secondAxis ? add_break!(data,BP[],pC) : BP[]
			bl1, dL   = support_fits!(data,deepcopy(pC), kdir, kL, breaks = b0)
			bline, dR = support_fits!(data,deepcopy(pC), kdir, kR,
			                          breaks= isempty(bl1)  ? b0 :  reverse!(bl1))
		end
		!iszero([dL,dR]) && break
		return t, bline
	end
	return 0,  BP[]
end



"""
Convolution at a particular point. Returns de filtered point
Parameters
==========
data: mesh & field data structure (FIELD)
eID : local element
zeta: relative position (in reference coordinates) at eID
type: either line or tensor

Optional
========
kx: 1st direction: list of number & order of bsplines for symmetric and shifted (default: 2p+1 splines of order p + 1)
ky: 2nd direction (only for tensor filter): list of number & order of bsplines for symmetric and shifted (default: 2p+1 splines of order p + 1)
theta: rotation in the kernel axis
scaling: kernel scaling (sized, ie, scaling = 1 means size 1, not size 1 * element_size )
"""
function filter_convolution(msh::MESH, fld::FIELD, eID::Int, zeta::vecD,type::String;
	 						kx=zeros(Int,0)::vecI, ky=zeros(Int,0)::vecI, theta = "auto"::Union{String,dbl,vecD},
							scaling = "auto"::Uniont{String,dbl,vecD}, inexact = false::Bool, track_footprint=false::Bool)
	p = fld.basis.degree
	if type == "tensor"
		isempty(kx) && (kx = [2*p+1,p+1,2*p+1,p+1,1])
		isempty(ky) && (ky = [2*p+1,p+1,2*p+1,p+1,1])
		if fld.basis.collapse
			if !msh.structured
				@error " Tensor product only supported for structured triangles"
				return 0.0
			end
			@error " Only periodic boundary conditions are considered. If the domain is non-periodic, near bounds will return unfiltered values"
			return tensor_tri_convolution(msh,fld,[KERNEL(kx[1],kx[2]),KERNEL(ky[1],ky[2])],eID, zeta, scaling)
		else
			return tensor_convolution(msh,fld,[varKernel(kx) varKernel(ky)],eID, zeta, theta, scaling, inexact,track_footprint)
		end
	else
		isempty(kx) && (kx = [2+p+1,p+1,2+p+1,p+1,1])
		return line_convolution(msh,fld,varKernel(kx),eID, zeta, theta,scaling, track_footprint)
	end
end


"""
This is a high-level function: computes the convolution for the tensor filter & quad elements
"""
function tensor_convolution(msh::MESH,fld::FIELD, kernels::Matrix{KERNEL}, eID::Int,
	                        zeta::vecD, theta::Union{vecD,String}, scaling::Union{vecD,String}, inexact::Bool, track_footprint::Bool)

	# FRAME : USE BDRY 1 AS "SOUTH" AND BDRY 2 AS "WEST"
	z1 = fld.zwC[1].nodes ; w1 = fld.zwC[1].weights
	z2 = fld.zwC[2].nodes ; w2 = fld.zwC[2].weights

	if theta == "auto" || scaling == "auto"
		vs = getfield.(msh.nodes[msh.elmts[eID].nodes],:xyz)
		he = msh.elmts[eID].h
		ca = CircularArray([1:msh.elmts[eID].type;])
		va = [acos(clamp(dot(vs[ca[j+1]].-vs[ca[j]],vs[ca[j]].-vs[ca[j-1]])/(he[ca[j+1]] * he[ca[j]]), -1.0, 1.0)) for j = 1:3]
		if msh.structured
			ra = findfirst(x->isapprox(x, 0.5 * pi), va)
			hx = abs(ndot(msh.frame[1,:],vs[ca[ra-1]].-vs[ca[ra]]))
			scaling = isapprox(hx,1.0) ? [he[ra-1],he[ca[ra]]] : [he[ra],he[ca[ra-1]]]
		else
			scaling = [maximum(he),maximum(he)]
		end
		kDIR = [msh.frame[1,:] .* scaling[1] msh.frame[2,:] .* scaling[2]]
	else
		kDIR = [[cos(theta[1]),sin(theta[1])] .* scaling[1] [cos(theta[2]),sin(theta[2])] .* scaling[2]]
	end
	centre = BP(standard_coords(zeta, element_coords(msh,eID,2)),[0.0,0.0], [eID],0)
	# MAIN AXIS
	kT = [eSYM, eSYM]
	kT[1], mA = get_element_list!(kernels[:,1], msh, centre, kDIR[1,:], skip_knots=inexact)
	isempty(mA) && (return evaluate_point(fld,eID,zeta))

	convolution = 0.0
	ix,iy = tensor_index(length(z1),length(z2))

	if track_footprint
		fig = figure(1)
		af = fig.gca()
		for j = 1 : msh.N
			ca = ele_ca(msh.elmts[j].type)
			coords = element_coords(msh,j)
			af.plot(coords[1,ca], coords[2,ca], c = "indigo")
			cen = element_centre(msh,j)
			#af.text(cen[1], cen[2], string(j))
		end
		fig2 = figure(2)
		ax = [fig2.add_subplot(121), fig2.add_subplot(122)]
	end

	for l = 1:length(mA) - 1
		aux = intersect(mA[l].e, mA[l+1].e)
		isempty(aux) && continue
		kT[2], b1 = get_element_list!(kernels[:,2], msh, BP(mA[l].z, mA[l].t, [aux[1]], mA[l].v), kDIR[2,:], secondAxis = true , skip_knots=inexact)
		isempty(b1) && (return evaluate_point(fld,eID,zeta))

		k2, b2 = get_element_list!(kernels[:,2], msh, BP(mA[l+1].z, mA[l+1].t, [aux[1]], mA[l+1].v), kDIR[2,:],kT = [kT[2]], secondAxis = true, skip_knots=inexact)
		isempty(b2) && (return evaluate_point(fld,eID,zeta))
		[kT[d] == ePOS && (kernels[ePOS,d].c = kernel_coefficients(kernels[ePOS,d].T, kernels[ePOS,d].zw)) for d = 1:2]
		lb1  = length(b1) ; lb2 = length(b2) ;
		c1   = 1 ; c2 = 1;
		supp_box =kernel_width(kernels[kT[2],2]) * scaling[2]
		while (c1 < lb1 || c2 < lb2)
			d1 = min(lb1, c1+1) ; 	d2 = min(lb2, c2+1)
			(d1 == c1 && d2 == c2 ) && break
			norm(b1[d1].z .- b1[c1].z) > supp_box && (d1 = c1)
			norm(b2[d2].z .- b2[c2].z) > supp_box && (d2 = c2)
			if d1 == c1 && d2 == c2
				c1 = min(lb1, c1+1) ; 	c2 = min(lb2, c2+1)
				continue
			end
			## for break points that are vertices, collect all surrounding elements
			ec1 = b1[c1].v >0 ? msh.nodes[b1[c1].v].ele : b1[c1].e
			ed1 = b1[d1].v >0 ? msh.nodes[b1[d1].v].ele : b1[d1].e
			ec2 = b2[c2].v >0 ? msh.nodes[b2[c2].v].ele : b2[c2].e
			ed2 = b2[d2].v >0 ? msh.nodes[b2[d2].v].ele : b2[d2].e

			## try finding common elements to see if we can integrate directly
			i1  = intersect(ec1, ed1) ;	i2  = intersect(ec2, ed2)
			i12 = intersect(i1,i2)
			if isempty(i12) && d1 > c1 && d2 > c2
				if !isempty(intersect(i1, ec2))
					d2 = c2 ;i2 = ec2
				elseif !isempty(intersect(i2, ec1))
					d1 = c1; i1 = ec1
				end
				i12 = intersect(i1,i2)
			end
			#draw kernel box in physical & kernel space
			kBOX = [b1[c1].t b2[c2].t]
			kXYZ = [b1[c1].z b2[c2].z]
			if d2 == c2 + 1
				kBOX = hcat(kBOX,b2[d2].t)
				kXYZ = hcat(kXYZ,b2[d2].z)
			end
			if d1 == c1 + 1
				kBOX = hcat(kBOX,b1[d1].t)
				kXYZ = hcat(kXYZ,b1[d1].z)
			end
			h12 = BP[]
			## if no intersection, draw a box to find all possible extra cuts with the mesh
			if isempty(i12)
				if isempty(intersect(ec1,ec2))
					p1  = BP(b1[c1].z, b1[c1].t, [b1[c1].e[1]], b1[c1].v)
					h12,dumb = support_fits!(msh, p1,b2[c2].z .- p1.z,[0.0,1.0])
				end
				if isempty(intersect(ed1,ed2))
					p2  = BP(b2[d2].z, b2[d2].t, [b2[d2].e[1]], b2[d2].v)
					h12,dumb = support_fits!(msh, p2, b1[d1].z .- p2.z,[0.0,1.0],breaks = h12)
				end
			end
			## Now add the bounding box as breaks for the integrals
			h12  = add_break!(msh,h12, b1[c1]) ; h12  = add_break!(msh,h12, b1[d1]) ;
			h12  = add_break!(msh,h12, b2[c2]) ; h12  = add_break!(msh,h12, b2[d2]) ;

			vl = findall(getfield.(h12,:v) .> 0)
			[h12[j].e = msh.nodes[h12[j].v].ele for j in vl ]
			totE = unique!(vcat(getfield.(h12,:e)...))
			for e in totE
				# look around the element vertices and check if vertex lives in kernel box
				for k in msh.elmts[e].nodes
					aux = findfirst(x->x == k, getfield.(h12,:v))
					if aux != nothing
						e ∉ h12[aux].e && push!(h12[aux].e,e)
						continue
					end
					zt = dbl[]
					if ccwV(kXYZ[:,1],kXYZ[:,2],kXYZ[:,3]) > 0 &&
						in_polygon(msh.nodes[k].xyz[1:2], kXYZ[:,1:3]) != -1
						z12 = local_coords(msh.nodes[k].xyz[1:2], kXYZ[:,1:3])
						zt  = standard_coords(z12, kBOX[:,1:3])
					elseif size(kBOX,2) == 4 && ccwV(kXYZ[:,1],kXYZ[:,3],kXYZ[:,4]) > 0 &&
						ccwV(kXYZ[:,1],kXYZ[:,3],kXYZ[:,4]) > 0 &&in_polygon(msh.nodes[k].xyz[1:2], kXYZ[:,[1,3,4]]) != -1
						z12 = local_coords(msh.nodes[k].xyz[1:2], kXYZ[:,[1,3,4]])
						zt = standard_coords(z12, kBOX[:,[1,3,4]])
					end
					!isempty(zt) && (h12 = vcat(h12,BP(msh.nodes[k].xyz[1:2], zt, [e], k)))
				end
				# Now collect all breaks corresponding to same element
				eles = findall(x-> e in x, getfield.(h12,:e))
				length(eles) < 3 && continue # not enough to integrate. SKIP !
				vs  = getfield.(h12[eles],:z)
				vID = convex_polygon(vs)
				ni  = size(vs,1)
				ca = CircularArray([1:ni;])
				####
				if track_footprint
					poly = hcat(vs[vID]...)
					af.fill(poly[1,ca[1:ni+1]], poly[2,ca[1:ni+1]], lw = 3, zorder = 0)
					af.plot(poly[1,ca[1:ni+1]], poly[2,ca[1:ni+1]], lw = 1, c = "red",zorder = 4)
				end
				###
				vO = vs[vID]
				n0 = 1
				normP   = [norm(vO[ca[j+1]] .- vO[ca[j]]) for j = 1:ni]
				pAngles = [dot(vO[ca[j-1]] .- vO[ca[j]],vO[ca[j+1]] .- vO[ca[j]])/(normP[ca[j-1]]*normP[ca[j]]) for j = 1:ni]
				pC = argmin(pAngles)
				while (n0 <= ni-2)  # triangulate polygon
					pv = ni == 4 && all( x-> abs(x.- pAngles[1]) .< fTOL, pAngles[2:end]) ? [1:4;] :
					                vcat(pC, ca[[pC+n0, pC+n0+1]])
					nv  = length(pv)
					sR  = hcat([local_coords(vO[j],element_coords(msh,e)) for j in pv]...)
					kR  = hcat([h12[eles[vID[j]]].t ./ scaling for j in pv]...)
					kXY,jK = standard_region(z1,z2,kR)
					kval = [evaluate_kernel(kernels[kT[j],j].T,kernels[kT[j],j].c, kXY[j,:]) for j = 1:2]
					if track_footprint
						cnv = ele_ca(nv)
						# kernel coordinates
						af.fill(poly[1,pv[cnv]], poly[2,pv[cnv]])
						ax[1].scatter(kXY[1,:], kval[1], ls = "--")
						ax[2].scatter(kXY[2,:], kval[2], ls = "--")
					end
					## get local coordinates (elements can be split)
					sXY,jac = standard_region(z1,z2,sR)  # coordinates for evaluating function
					pLeg = eval_basis(fld.basis,sXY[1,:], sXY[2,:])
					proj = [sum(fld.modes[:, e] .* pLeg[i,:]) for i = 1:length(ix) ]
					convolution += sum(proj .* kval[1] .* kval[2] .* w1[ix] .* w2[iy] .* jK)
					n0 += nv - 2
				end
				ni == length(h12) && break # already done loop around any more breaks in h. all have been accounted for
			end
			c1  = d1; c2 = d2
		end
	end
	return convolution
end

"""
This is a high-level function: computes the convolution for the line filter for quad & tris elements
"""
function line_convolution(msh::MESH,fld::FIELD, kernel::Vector{KERNEL}, eID::Int, zeta::vecD,
	                      theta::Union{dbl,String},scaling::Union{dbl,String}, track_footprint::Bool)

	z   = fld.zwC[1].nodes ; w = fld.zwC[1].weights
	# MAIN AXIS
    if theta == "auto" || scaling == "auto"
        if mesh_type(msh) == 3
            he = msh.elmts[eID].h
			ca = ele_ca(3)
			vj = [msh.nodes[msh.elmts[eID].nodes[ca[j+1]]].xyz .- msh.nodes[msh.elmts[eID].nodes[ca[j]]].xyz for j = 1:3]
			va = [acos(clamp( dot(vj[ca[j+1]],vj[ca[j]])/(he[ca[j+1]] * he[ca[j]]), -1.0, 1.0)) for j = 1:3]
            (theta == "auto") && (theta = minimum(vcat(va, pi .- va)))
			(scaling == "auto") && (scaling = maximum(he))

        else
            h1 = msh.elmts[eID].h[1]
            h2 = msh.elmts[eID].h[2]
            if scaling == "auto"
                scaling = theta == 0.0 ? h1 : theta == 0.5 * pi ? h2 : sqrt(h1^2+h2^2)
			end
			theta == "auto" && (theta = atan(h2/h1))
        end
    end

	kDIR   = [cos(theta), sin(theta)] .* scaling
	centre = BP(standard_coords(zeta, element_coords(msh,eID,2)),[0.0,0.0], [eID],0)
	kT,kB  = get_element_list!(kernel, msh, centre, kDIR)

	if isempty(kB)
		theta = pi - theta
		kDIR = [cos(theta), sin(theta)] .* scaling
		kT,kB = get_element_list!(kernel, msh, centre, kDIR)
		if isempty(kB)
			theta = pi - theta
			scaling *= 0.5
			kDIR = [cos(theta), sin(theta)] .* scaling
			kT,kB = get_element_list!(kernel, msh, centre, kDIR)
			if isempty(kB)
				kDIR = [cos(theta), sin(theta)] .* scaling
				kT,kB = get_element_list!(kernel, msh, centre, kDIR)
				isempty(kB) && return evaluate_point(fld,eID,zeta)
			end
		end
	end

	kT == ePOS && (kernel[ePOS].c = kernel_coefficients(kernel[ePOS].T, kernel[ePOS].zw))
	convolution = 0.0

	if track_footprint
		fig = figure(3)
		af = fig.gca()
		for j = 1 : msh.N
			ca = ele_ca(msh.elmts[j].type)
			coords = element_coords(msh,j)
			af.plot(coords[1,ca], coords[2,ca], c = "indigo")
			cen = element_centre(msh,j)
			#af.text(cen[1], cen[2], string(j))
		end
		for j = 1 : length(kB)
			af.scatter(kB[j].z[1], kB[j].z[2], c = "r", zorder = 4)
			j == length(kB)-1 && break
			af.plot([kB[j+1].z[1],kB[j].z[1]],[kB[j+1].z[2],kB[j].z[2]])
		end
	end

	for l = 1:length(kB)-1
		### KERNEL COORDINATES t: breaks in kernel space ( t = b/H )
		le = intersect(kB[l].e,kB[l+1].e)
		isempty(le) && continue

		t2 = norm(kB[l+1].t) / scaling * sign(dot(kDIR,kB[l+1].t));
		t1 = norm(kB[l].t)   / scaling * sign(dot(kDIR,kB[l].t));

		dk = 0.5 * (t2-t1)#nt   # parametric space
		# Map for gauss quadrature
		kerT  = dk.* z .+ 0.5 * (t2 +t1)
		kval  = evaluate_kernel(kernel[kT].T,kernel[kT].c, kerT)
		### DG COORDINATES z (phys): map to local (x) -> chi -> eta
		lxy   = [0.5.*((kB[l+1].z[j] - kB[l].z[j]) .* z .+ kB[l+1].z[j] .+ kB[l].z[j]) for j = 1:2]
		ve    = element_coords(msh,le[1])
		legXY = hcat([local_coords([lxy[1][j],lxy[2][j]],ve) for j = 1:length(z)]...)
		pLeg  = eval_basis(fld.basis,legXY[1,:], legXY[2,:])
		proj  = [sum(fld.modes[:, le[1]] .* pLeg[i,:]) for i = 1:length(z)]
		convolution  += sum(proj .* kval .* w) * dk
	end

	return  convolution
end

"""
Computes the L2 errors for the entire field
Parameters:
===========
h: reference size (scalar ( h1 * h2))
f1: first field
f2: second field
w1,w2: Gauss weights in each direction
"""
function get_eL2(h::dbl, f1::matD, f2::matD, w1::vecD, w2::vecD)
	q,n   = size(f1)
	i1,i2 = tensor_index(length(w1), length(w2))
	if q != length(i1)
		@info "get_eL2 !!! miss-match between data and quadrature weigths!! f1 has $q points (nodes) and we have $(length(i1)) weights"
		return 0.0
	end
	return sqrt(0.25 * h * sum([sum((f1[:,e] .- f2[:,e]).^2 .* w1[i1] .* w2[i2]) for e = 1:n]))
end


"""
Test function: filters a particular data point
Parameters:
===========
rx: number of splines in the 1st direction / line direction
lx: order of splines in the 1st direction / line direction
ry: number of splines in the 2nd direction (for line filter, write 0)
ly: order of splines in the 2nd direction (for line filter, write 0)
data: mesh & field data structure (twoD_mesh_field)
e: element id
zeta: reference coordinates at element
kType: kernel type (line / tensor)

Optional:
========
theta = kernel orientation (default Cartesian for tensor and diagonal for line)
scaling = kernel scaling (default element size (hx,hy) and element diagonal (for line))
"""
test_filter(rx::Int, lx::Int, ry::Int, ly::Int, msh::MESH,fld::FIELD, e::Int,zeta::vecD,kType::String;
	                 theta = "auto"::Union{vecD,dbl,String}, scaling="auto"::Union{vecD,dbl,String}, track_footprint=false::Bool) =
                     filter_convolution(msh,fld, e, zeta , kType, track_footprint=track_footprint,
					 kx =[rx,lx,rx,lx,1], ky =[ry,ly,ry,ly,1], theta = theta, scaling=scaling)

test_filter(rx::Int, lx::Int, ry::Int, ly::Int, data::DATAFIELD, e::Int,zeta::vecD,kType::String;
	                 theta = "auto"::Union{vecD,dbl,String}, scaling="auto"::Union{vecD,dbl,String}, track_footprint=false::Bool) =
                  filter_convolution(data.mesh,data.fields[data.fdict[""]], e, zeta, kType, track_footprint=track_footprint,
					 kx =[rx,lx,rx,lx,1], ky =[ry,ly,ry,ly,1], theta = theta, scaling=scaling)

end

