
#=
*
*      MSIAC: the SIAC kernel block.
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
=#
const eSYM  = 1
const ePOS  = 2
const BS = BSplines

export plot_kernel, test_kernel

"""
SIAC KERNEL K(r+1,l). Members:
-----------------------
r = total splines -1

l = spline order

c = kernel coefficients vector (zeros for type ePOS)

type = eSYM (symmetric) OR ePOS (shifted)
"""
mutable struct KERNEL
	l  :: Int
	rC :: Int
	rB :: Int
    T  :: matD
	c  :: vecD
	zw :: GZW
	function KERNEL(rC::Int,l::Int ; rB=0::Int, type=eSYM::Int)
		T  = knot_matrix(l, rC; rB = 0, lX = 0.0)
		l == 1 && (rB = 0 )
		zw = GZW("legendre", l + rC + rB)
        c  = type == eSYM ? kernel_coefficients(T,zw) : zeros(rC + rB)
		return new(l, rC, rB, T, c, zw)
    end
end

# rs ls conditions for symmetric rp lp conditions for shifted kernel
varKernel(pars::vecI) = [KERNEL(pars[1],pars[2]),KERNEL(pars[3],pars[4]; rB = pars[5], type = ePOS)]
kernel_width(k::KERNEL) = maximum(k.T) - minimum(k.T)
function kernel_coefficients(T::matD, zw::GZW)
	m,n = size(T)
	# Compute the kernel coefficients
	A  = zeros(m,m)
	# SHIFTED COEFFICIENTS X #
	z0 = 0.0
	yp = vcat(1.0,zeros(m-1))
	zQ = zw.nodes ; wQ = zw.weights
	for p = 1:m
		for g = 1:m
			spline = BS.BSplineBasis(n-1, T[g,:])[n-1]
			for k = 1:n-1
				baM = T[g,k] - T[g,k+1] ; baP = T[g,k] + T[g,k+1]
				y   = 0.5 .* ( baM .* zQ .+ baP)
				A[p,g] -= sum((z0.-y).^(p-1) .* BS.splinevalue.(spline,y) .* wQ) * baM * 0.5
			end
		end
	end
	c  = factorize(A)\yp
	#!(isapprox(sum(c),1.0,atol=1.e-14)) && @warn " C  ",c,"  SUM (C_gamma) =$(sum(c)) != 1  ERR = $(sum(c)-1) !!"
	return c
end

"""
Returs the c_gammas. If lambda != 0, computes the shifted coeffs.
Otherwise, returns stored symmetric coefficients.
"""
function analytic_kernel_coeff(r::Int, l::Int, lambda=0.0::dbl)
    A  = zeros(r+1,r+1)
    # SHIFTED COEFFICIENTS X #
    for m = 1:r+1
        s1 = [sum([(-1.0)^(j+l-1) * dbl(binomial(l-1,j)) *
                ((dbl(j)-0.5*(dbl(l) -2.0) )^(l+n-1) - (dbl(j)-0.5*dbl(l))^(l+n-1)) for j = 0:l-1]) for n = 1:m]
        for gam = 1:r+1
            xShift   = dbl(gam - 1) + lambda - 0.5 * r
            s2       = [binomial(m-1,n-1)*(xShift)^(m-n)*dbl(factorial(n-1)/factorial(n-1+l)) for n = 1:m]
            A[m,gam] = sum(s1 .* s2)
        end
    end
    #return reverse(factorize(A)\vcat(1.0,zeros(r)))
	return factorize(A)\vcat(1.0,zeros(r))
end

"""
Returns the knot matrix T based on r, l and shfit
"""
function knot_matrix(l::Int, rC::Int; rB =0::Int, lX=0.0::dbl)
	lK = lX < 0.0 ? lX + (rC+l-1)*0.5 : lX - (rC+l-1)*0.5
	T  = [ j+i-0.5*(rC+l-1)+lX for i = 0:rC-1, j = 0:l]
	#@info " SPLINES C ",rC, " B ",rB," FULL SHIFT ", lX, " (x-xB)/H ", lK
	if rB > 0 && lX != 0.0
		Tg = lX < 0.0 ? vcat(lK-1, fill(lK, l) )' : vcat(fill(lK,l), lK+1)'
		T  = lX < 0.0 ? vcat(T, Tg) : vcat(Tg, T)
	end
	#@info " KNOT MATRIX ",T
	return T
end

"""
Update knot matrix: useful when needed reshifting
"""
function update_knot_matrix!(ker::KERNEL ; lX=0.0::dbl)
	ker.T =  knot_matrix(ker.l, ker.rC; rB = ker.rB, lX=lX)
end


"""
Evaluate kernel at given points
Params
------
ker:
    KERNEL
z:
    sample points
"""
function evaluate_kernel(T::matD, c::vecD, z::vecD)

    #[z[i]==-0.0 && (z[i] = 0.0) for i=1:length(z)]
	m,n = size(T)
    # Obtain the splines value at the gauss points
    splines = zeros(m,length(z))
    for j = 1:m
        local bs = BS.BSplineBasis(n-1, T[j,:])[n-1]
        splines[j,:] = BS.splinevalue.(bs, z)
    end
    return [sum(splines[:,k] .* c) for k = 1:length(z)]
end


"""
Given a knot matrix, sorts it and splits it into >0 and <0
Params
------
T:
    knot matrix
"""
function split_knots(T::matD)
	uT = unique!(vcat(T...))
	nI = filter(x -> x < 0.0, uT)
	pI = filter(x -> x >= 0.0, uT)
	return reverse(sort(nI)), sort(pI)
end


"""
Plots kernel and spline functions
Params
------
c:
    kernel coefficients
T:
    knot matrix
"""
function plot_kernel(T::matD, c::vecD)
    k_supp = [minimum(T), maximum(T)]
    N      = 1000
    x      = [k_supp[1]:(k_supp[2] - k_supp[1]) / N : k_supp[2];]
    m,n = size(T)
    # Obtain the splines value at the gauss points
    splines = zeros(m,length(x))
    figure(1, figsize=(2,2))
    c_s = ["red", "orange"]
	l_s = ["--","-."]
    c_k = "black"
    for j = 1:m
        local bs     = BS.BSplineBasis(n-1, T[j,:])
        splines[j,:] = BS.splinevalue.(bs[n-1], x)
        piv = findall(x -> !iszero(x), splines[j,:])
		lab = j == 1 ? "BS_"*string(n-1) : ""
        plot(x[piv], splines[j,piv], c = c_s[j % 2 + 1], ls = l_s[j%2 + 1], label = lab)
    end
    ker  = [sum(splines[:,k] .* c) for k = 1:length(x)]
    plot(x, ker, c = c_k, ls = "-", label = "K("*string(m)*","*string(n-1)*")")
    legend()
    show(block=true)
end


function test_kernel(k = 2::Int; shift = 0.0, genSpline = 1)
    r   = 2*k+1 ; l = k+1
	ker = KERNEL(r, l, rB=genSpline)
	lambda = -0.5*shift*(3*k+1)
	update_knot_matrix!(ker, lX=lambda)
	ker.c = kernel_coefficients(ker.T, ker.zw)
	plot_kernel(ker.T,ker.c)
end
