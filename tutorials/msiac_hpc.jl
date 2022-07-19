
using DelimitedFiles, Printf, Distributed
include(ENV["MSIAC"]*"/src/MSIAC_HPC.jl")


function test_distr(sim, p, ft)
	nE        = 40
	domX      = [0.0,2.0]
	domY      = [0.0,2.0]
	solType   = "sine"
	CPUtimes  = 0.0
	np        = nworkers()
	println(" USING THREADS $np  total points $(nE^2 * 25)")
        mf = create_2D_mesh(4,nE,nE,domX,domY,true,"mymesh.txt")
	ff = create_analytic_field([p+1,p+1],["legendre","legendre"],mf,solType*"_fld.txt",type=solType)
	data = load_data(mf,ff); 
	update_phys!(data, "Legendre",5,"Legendre",5)
	dg     = data.fields[1].f
	p1, p2 = get_evaluation_points(data)
	w1, w2 = get_evaluation_weights(data)
	exact  = get_exact_data(data.mesh,solType,p1,p2)  #we are cheating a bit...
	dgL2   = get_eL2(characteristic_size(data.mesh), exact, dg,w1,w2)

	CPUtimes = @elapsed filter = filter_data_hpc(data,ft)
	siacL2  = get_eL2(characteristic_size(data.mesh), exact, filter,w1,w2)

	println("\n\n ========= RESULTS FOR p = $p ================= \n")
	println("N\t ||dg e||_0\t\t  ||Distr||\t\t\n")
	println("------------------------------------------------------\n")
	@printf("%3s x%3s\t %.2e\t %.2e\t %lf\t \n",
               nE,nE, dgL2,siacL2, CPUtimes)

	println("------------------------------------------------------\n")
	return CPUtimes
end


for kDG = 1 : 3
	t1 = @elapsed f1 = test_distr(1,kDG, "line")
	@info " RUNNING FOR $(nprocs()) k DG $kDG  TOTAL TIME $t1 FILTERING $f1 "
end
