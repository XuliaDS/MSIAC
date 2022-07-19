#=
*
*      MSIAC: unit testing
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
=#
push!(LOAD_PATH, ".")
using Test, Printf, MSIAC
locFdr = joinpath(dirname(pathof(MSIAC)),"../test/")
### Generate mesh & field files ####
function generate_meshes()
	@info " Mesh files not found. Generating 5 files "
	aux = create_2D_mesh(4,7,7,  [0.0,2.0],[0.0,2.0],  false, locFdr*"mesh1.txt")
	aux = create_2D_mesh(4,12,12,[0.0,2.0],[0.0,2.0],  false, locFdr*"mesh2.txt")
	aux = create_2D_mesh(4,12,12,[0.0,1.0],[0.0,1.0],  false, locFdr*"mesh3.txt",pert = [0.25,0.4])
	aux = create_2D_mesh(3,8,8,  [0.0,2.0],[0.0, 2.0], true , locFdr*"mesh4.txt")
	aux = create_2D_mesh(3,8,8,  [0.0,2.0],[0.0, 2.0], false, locFdr*"mesh5.txt",pert = [0.3,0.35], structured = false)
	@info " Done generating mesh files "
end
function generate_fields()
	@info " Field files not found. Generating 5 files "
	lg = "legendre" ; lb = "lobatto"
	aux = create_analytic_field([2,2],[lg,lg],locFdr*"mesh1.txt",locFdr*"field1.txt",type="poly_xy")
	aux = create_analytic_field([2,2],[lg,lg],locFdr*"mesh2.txt",locFdr*"field2.txt",type="poly_xy2")
	aux = create_analytic_field([2,2],[lb,lg],locFdr*"mesh3.txt",locFdr*"field3.txt",type="poly_xy2")
	aux = create_analytic_field([1,1],[lg,lg],locFdr*"mesh4.txt",locFdr*"field4.txt",type="ctt")
	aux = create_analytic_field([3,3],[lg,lg],locFdr*"mesh5.txt",locFdr*"field5.txt",type="poly_x+y")
	@info " Done generating field files "
end
# load data from mesh & field
load_mf_data(n::Int, exp::String, structured=false) = load_data(locFdr*"mesh"*string(n)*".txt",locFdr*"field"*string(n)*".txt"; modalExp=exp, structured=structured)

function test1_regular_quads_poly_reproduction_legendre()
	@info " \n\n Regular quads polynomial reproduction: x * y with legendre expansions "
    data  = load_mf_data(1,"legendre")
	p     = polynomial_degree(data)
	Q1    = 2 * p+1 ; Q2 = 2 * p + 2
	update_phys!(data, "lobatto",Q1,"legendre",Q2)
	fL    = filter_data(data,"line")
	fT    = filter_data(data,"tensor")
	fA    = get_field(data)
	eL    = maximum(abs.(fA .- fL)) ; eT = maximum(abs.(fA .- fT)) ;
    @info " MAX ERROR APPROXIMATION line $eL tensor $eT"
    return eL < 1.e-12 && eT < 1.e-12
end

function test2_regular_quads_poly_reproduction_hierarchy()
    @info " \n\n Regular quads polynomial reproduction: x * y^2 with filter and hierarchical expansion "
    data  = load_mf_data(2,"hierarchy")
	p1,p2 = get_evaluation_points(data)
	exact = get_exact_data(data, "poly_xy3", p1,p2)
	eD    = maximum(abs.(exact .-get_field(data)))
	fL    = filter_data(data,"line")
	eL    = maximum(abs.(exact .- fL))
	fT    = filter_data(data,"tensor")
	eT    = maximum(abs.(exact .- fT))
    @info " MAX ERROR APPROXIMATION line $eL tensor $eT  data $eD"
    return eL < 1.e-10 && eT < 1.e-10
end

function test3_irregular_quads_poly_reproduction_legendre()
    @info " \n\n Irregular quads polynomial reproduction: x * y^2 with Legendre expansion "
    data = load_mf_data(3,"legendre")
	p = polynomial_degree(data)
    Q1   = 2 * p -1 ; Q2 = 2 * p
	update_phys!(data, "lobatto",Q1,"legendre",Q2)
    fL = filter_data(data,"line")
    fT = filter_data(data,"tensor")
	fA = get_field(data)
	eL = maximum(abs.(fA .- fL)) ; eT = maximum(abs.(fA .- fT))
    @info " MAX ERROR APPROXIMATION line $eL  tensor $eT"
    return eL < 1.e-12 && eT < 1.e-12
end

function test4_line_structured_tris_poly_reproduction()
	@info " \n\n Structured triangles polynomial reproduction: ctt with Pk expansion "
	data = load_mf_data(4,"Pk",true)
	update_phys!(data, "Legendre",3,"legendre",3)
	fA   = get_field(data)
	fT   = filter_data(data, "tensor")
	fL   = filter_data(data, "line")
	eT   = maximum(abs.(fA .- fT)) ; eL = maximum(abs.(fA .- fL))
	@info " MAX ERROR APPROXIMATION line $eL tensor $eT "
	return eL < 1.e-12 && eT < 1.e-12
end

function test5_line_unstructured_tris_poly_reproduction()
    @info " \n\n Unstructured triangles Polynomial reproduction: x * y with Pk & Hierarchy expansions and line filters"
    dataP = load_mf_data(5, "Pk")
	dataH = load_mf_data(5, "hierarchy")
	p = polynomial_degree(dataP)
	Q1    = 2 * p +4 ; Q2 = Q1
	tp    = ["legendre","legendre"]
	update_phys!(dataP, tp[1],Q1,tp[2],Q2)
	update_phys!(dataH, tp[1],Q1,tp[2],Q2)
    fP = filter_data(dataP,"line")
	fH = filter_data(dataH,"line")
	eP = maximum(abs.(get_field(dataP) .- fP)) ;eH = maximum(abs.(get_field(dataH) .- fH)) ;
    @info " MAX ERROR APPROXIMATION Pk $eP  Hierarchy $eH  "
	return eP < 1.e-12 &&  eP < 1.e-12
end

function test_convergence_sine_field(nEX,nEY, nR, pmin, pmax; type = "sine")
	@info "\n\n Superconvergence test for a sine wave over a uniform quad mesh"
	domX = [0.0, 2.0]
	domY = [0.0, 2.0]
	for o = 0:pmax-pmin
        p     = pmin + o
        eL2   = zeros(nR + 1, 3)
        eLInf = zeros(nR + 1, 3)
        eleX  = zeros(Int64,nR + 1)
		eleY  = zeros(Int64,nR + 1)
		Q     = p + 2
	    for k = 1:nR+1
            eleX[k] = 2 ^(k - 1) * nEX
			eleY[k] = 2 ^(k - 1) * nEY
			# create mesh and field :: periodic BCs
			myMesh = create_2D_mesh(4,eleX[k],eleY[k],domX,domY,false, "mymesh.txt";
			                        yFirst = false, pert = [0.0,0.0], structured = true)
			field  = create_analytic_field([p+1,p+1], ["legendre","legendre"], myMesh, type*"_field.txt";type=type)
			data   = load_data(myMesh, field, modalExp = "legendre"; degree = p);
			data.zwS[1] = gaussZW("legendre",Q); data.zwS[2] = gaussZW("legendre",Q)
			# FILTER FIELD: WILL REPROJECT ONTO A MODAL p BASIS
			 fL = filter_data(data,"line")
			 fT = filter_data(data,"tensor")
			# CREATE A FAKE MODAL SOLUTION: PROJECT FIELD ONTO P BASIS
			approx = update_field(data)
			# LOAD DATA FOR PP: error analysis
			exact = reshape(eval_analytic_field(mesh_local_to_global(data);type=type),data.zwS[1].Q*data.zwS[2].Q,data.N)  # transform into element matrix
			# GET L2 and LINF ERRORS: CONVERGENCE ANALYSIS
			h = hcat(getfield.(data.elmts,:h)...)
			eL2[k,:]   = [sqrt(get_eL2(h, exact,approx,data.zwS[1].weights,data.zwS[2].weights))  ,#/ (sX*sY),
			              sqrt(get_eL2(h, exact,fL    ,data.zwS[1].weights,data.zwS[2].weights)) ,
						  sqrt(get_eL2(h, exact,fT    ,data.zwS[1].weights,data.zwS[2].weights))]#/ (sX*sY)]
			eLInf[k,:] = [maximum(abs.(exact .- approx)),maximum(abs.(exact .- fL)),maximum(abs.(exact .- fT))]
			rm(myMesh) ; rm(field) ;
        end
    	@info "\n\n =========================== RESULTS FOR p = $p =============================== "
        @info "N\t ||e||_0\t ORDER\t ||e_Line*||_0\t ORDER\t ||e_Tensor*||_0\t ORDER\t ||e||_inf\t ORDER\t ||e_Line*||_inf\t ORDER \t||e_Line*||_inf\t ORDER"
        @info "-----------------------------------------------------------------------------------"
        for k in 1:nR+1
			ne = string(eleX[k])*" x "*string(eleY[k])
            if k == 1
                @printf("%9s\t %.2e\t NA \t %.2e\t NA  \t %.2e\t NA  \t %.2e\t NA \t %.2e\t NA \t %.2e\t NA\n",ne, eL2[k,1],eL2[k,2],eL2[k,3],eLInf[k,1],eLInf[k,2],eLInf[k,3])
            else
                r2 = log10.(eL2[k-1,:]./eL2[k,:])./log10(2)
                ri = log10.(eLInf[k-1,:]./eLInf[k,:])./log10(2)
				@printf("%9s\t%.2e\t %.2f\t %.2e\t %.2f\t %.2e\t %.2f\t %.2e\t %.2f\t %.2e\t %.2f %.2e\t %.2f\t\t\n",
                         ne, eL2[k,1],r2[1], eL2[k,2],r2[2], eL2[k,3],r2[3], eLInf[k,1],ri[1],eLInf[k,2],ri[2],eLInf[k,3],ri[3] )
            end
        end
		print("-----------------------------------------------------------------------------------\n")
	end
	return true
end

@testset " SIAC2D API Unit Testing " begin
	nt = 5
	for j = 1: nt
		mf = "mesh"*string(j)*".txt"
		ff = "field"*string(j)*".txt"
		if !isfile(locFdr*mf) || !isfile(locFdr*ff)
			generate_meshes()
			generate_fields()
			break
		end
	end
	@test test1_regular_quads_poly_reproduction_legendre()
	@test test3_irregular_quads_poly_reproduction_legendre()
    @test test4_line_structured_tris_poly_reproduction()
    @test test5_line_unstructured_tris_poly_reproduction()
	#@info " Higher order &  Superconvergence test ? (this will take > 10 minutes ) 1 = yes, 0 = no"
	#if readline() == "1"
	#	@test test2_regular_quads_poly_reproduction_hierarchy()
	#	@test test_convergence_sine_field(10,10,2,2,3)
	#end
end
