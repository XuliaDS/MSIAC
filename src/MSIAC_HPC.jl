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

using Distributed 
@everywhere push!(LOAD_PATH,ENV["MSIAC"])
 using MSIAC
@everywhere using Distributed, DistributedArrays 


export filter_data_hpc

@everywhere function parallel_conv!(res, eC, i1, i2, msh, fld,kType, kx,ky,theta,scaling,inexact)
	for e in 1:length(eC)
		res[e] = filter_convolution(msh,fld,eC[e],[i1[e],i2[e]],
		                            kType,kx=kx,ky=ky,theta=theta,scaling=scaling,inexact=inexact)
	end
end




function filter_data_hpc(msh::MESH, fld::FIELD, kType::String; rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
	                 dumpFile="nofile"::String, theta="auto"::Union{String,Float64,Vector{Float64}}, scaling="auto"::Union{String,Float64,Vector{Float64}},
					 inexact=false::Bool)

	(rx == -1) && (rx = 2 * fld.basis.degree + 1)
	(ry == -1) && (ry = 2 * fld.basis.degree + 1)
	(lx == -1) && (lx =     fld.basis.degree + 1)
	(ly == -1) && (ly =     fld.basis.degree + 1)

	Q = max(2,ceil(Int, 0.5 * (3 * fld.basis.degree + max(lx, ly) + 1)))
	fld.zwC[1] = GZW("legendre",Q)
	fld.zwC[2] = fld.basis.collapse ? GZW("Radau",Q) : GZW("legendre",Q)
	i1, i2  = tensor_index(fld.zw[1].Q,fld.zw[2].Q,fld.reverse)
	Qt   = fld.zw[1].Q * fld.zw[2].Q
	filt = zeros(Qt,msh.N)
	i1C  = vcat(repeat(fld.zw[1].nodes[i1], msh.N)...)
	i2C  = vcat(repeat(fld.zw[2].nodes[i2], msh.N)...)
	eC   = vcat([fill(e,Qt) for e in 1:msh.N]...)

	aux   = dzeros((Qt*msh.N,), workers(), nworkers())
	d_i1C = distribute(i1C, procs = workers())
	d_i2C = distribute(i2C, procs = workers())
	d_eC  = distribute(eC , procs = workers())
	@sync begin
		for p in workers()
			@spawnat p parallel_conv!(localpart(aux),localpart(d_eC),localpart(d_i1C),localpart(d_i2C),
				                          msh,fld,kType,[rx,lx,rx,lx,1],[ry,ly,ry,ly,1],theta,scaling,inexact)
		end
	end	
	filt = reshape(convert(Array, aux), Qt, msh.N)
	if dumpFile != "nofile"
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
	end
	#showWarnings && @info " Finished postprocessing field with $kType kernel. Total points: $(data.N * data.zw[1].Q * data.zw[2].Q) "
	@everywhere GC.gc()
	return filt
end
filter_data_hpc(data::DATAFIELD, kType::String; f=""::String,rx =(-1)::Int, lx=(-1)::Int, ry=(-1)::Int, ly=(-1)::Int,
	                 dumpFile="nofile"::String, theta="auto"::Union{String,Float64,Vector{Float64}},
					 scaling="auto"::Union{String,Float64,Vector{Float64}}, parallel=true::Bool, inexact=false::Bool) =
					 filter_data(data.mesh,data.fields[data.fdict[f]],kType, rx = rx, lx = lx, ry = ry, ly=ly,
					 dumpFile = dumpFile, theta=theta, scaling=scaling, parallel=parallel, inexact=inexact)


