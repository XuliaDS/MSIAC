#=
*
*      MSIAC: structured triangles. Convolution using tensor product kernels
*
*      Copyright 2022, Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Matthew Picklo
*      Email: mpicklo@mines.edu
*
=#


function tensor_tri_convolution(msh::MESH,fld::FIELD, kernels::Vector{KERNEL}, e, zeta, scaling)

	boundaries_x = hcat(getfield.(msh.nodes[msh.bMap[1]],:xyz)...)[1,:]
	boundaries_y = hcat(getfield.(msh.nodes[msh.bMap[3]],:xyz)...)[2,:]
	scaling == "auto" && (scaling=[boundaries_x[2] - boundaries_x[1],boundaries_y[2] - boundaries_y[1]])
	periodicity=msh.periodic
    Nx = length(boundaries_x) - 1
    Ny = length(boundaries_y) - 1

    order = fld.basis.degree + 1

    # set up criteria
    p = order - 1
    BSorderx = kernels[1].l
    RSx = kernels[1].rC
    BSordery = kernels[2].l
    RSy = kernels[2].rC


    Hx = scaling[1] ; Hy = scaling[2]
	hx = boundaries_x[2]   - boundaries_x[1] ; hy = boundaries_y[2]   - boundaries_y[1]
    Lx = boundaries_x[end] - boundaries_x[1]
    Ly = boundaries_y[end] - boundaries_y[1]

    #  NUMBER OF GAUSS POINTS FOR EVALUATION OF POLYNOMIAL
    #kgl = ceil(Int, 0.5 * (3 * p + max(BSorderx, BSordery) + 1))
    zl = fld.zwC[1].nodes
	wl = fld.zwC[1].weights   # need enough quad points to support k + BSOrder degree
	zr = fld.zwC[2].nodes
	wr = fld.zwC[2].weights   # need enough quad points to support k + BSOrder degree

    zx, zy = tensor_index(length(zl))
    wx, wy = tensor_index(length(zr))

    xi1, xi2 = collapsed_to_cartesian(zl[zx], zr[zy])


    #  Kernel stuff ##############################
    halfWidthX = Hx * 0.5 * (BSorderx + RSx - 1)
    halfWidthY = Hy * 0.5 * (BSordery + RSy - 1)


    kBrksTx = [-halfWidthX:Hx:halfWidthX;]
    kBrksTy = [-halfWidthY:Hy:halfWidthY;]

    NbrkTx = length(kBrksTx)
    NbrkTy = length(kBrksTy)

    R = zeros(4, 2, (NbrkTy - 1) * (NbrkTx - 1))
    #Need to break support up into rectangular regions defined by vertices in
    #CC order. These are in terms of the t-variables
    for j = 1:NbrkTy-1
        for i = 1:NbrkTx-1
            R[:, :, i+(j-1)*(NbrkTx-1)] = [kBrksTx[i] kBrksTy[j]
                kBrksTx[i+1] kBrksTy[j]
                kBrksTx[i+1] kBrksTy[j+1]
                kBrksTx[i] kBrksTy[j+1]]
        end
    end

    #Filterpoint location in global coords
    offset = standard_coords_stt(zeta, element_coords(data, e))

    #If periodic boundaries and not assummed and the kernel overlaps the boundaries, do not filter
    if !periodicity
        if (offset[1] - halfWidthX < boundaries_x[1]) || (offset[1] + halfWidthX > boundaries_x[end]) || (offset[2] - halfWidthY < boundaries_y[1]) || (offset[2] + halfWidthY > boundaries_y[end])
            return evaluate_point(data, e, zeta)
        end
    end

    #Determine super-element containing (x,y) (Only works over semi-uniform meshes)
    I, J = findQuadSuperElmnt(offset[1], offset[2], boundaries_x, boundaries_y)

    #Determine super-elements overlapped by kernel support
    if hx < Hx || hy < Hy
        overlapMat = kFootPrint(I, J, halfWidthX, halfWidthY, hx, hy, Nx, Ny)
    else
        overlapMat = kFootPrint2(halfWidthX, halfWidthY, boundaries_x, boundaries_y, offset[1], offset[2])
    end
    #Have to turn into integers, else Julia complains about the type
    overlapMat = Int.(overlapMat)
    overlapMat = unique(overlapMat, dims=1)

    R_coords = zeros(4, 2, (NbrkTy - 1) * (NbrkTx - 1))
    loc_array = []
    R_array = []
    count = 1
    for ii = 1:size(R, 3)

        #Define Kernel support in global coords
        pnt_filt = [offset[1] offset[2]]

        x1 = R[1, 1, ii]
        x2 = R[2, 1, ii]
        y1 = R[1, 2, ii]
        y2 = R[3, 2, ii]

        VNW = pnt_filt .+ [x1 y2]
        VSW = pnt_filt .+ [x1 y1]
        VNE = pnt_filt .+ [x2 y2]
        VSE = pnt_filt .+ [x2 y1]

        #Return vertices of separate kernel support polygons when support
        #overlaps boundaries. Note this supposes the kernel support is at
        #most as large as the domain itself.

        R_split = kspptPeriodic(boundaries_x[1], boundaries_x[end], boundaries_y[1], boundaries_y[end], VNW, VSW, VSE, VNE)

        amnt = size(R_split, 3)
        R_coords[:, :, ii] = [VNW; VSW; VSE; VNE]
        counter = 1
        for cntr = count:count+amnt-1
            loc_array = append!(loc_array, ii)
            if ii == 1
                R_array = R_split
            else
                R_array = cat([R_array, R_split[:, :, counter]]..., dims=3)
            end
            counter = counter + 1
        end

        count = count + amnt

    end
    #Loop over elements contained in overlapMat and evaluate integrals
    #where applicable
    intVal = 0

    for ii = 1:size(overlapMat, 1)
        #Determine lower and upper triangle vertices based off of I,J
        #coords

        I = overlapMat[ii, 1]
        J = overlapMat[ii, 2]

        KL, KR = superElmtToTri(I, J, Nx)

        L_nodes = msh.elmts[KL].nodes
        tL = hcat(getfield.(msh.nodes[L_nodes], :xyz)...)'

        R_nodes = msh.elmts[KR].nodes
        tR = hcat(getfield.(msh.nodes[R_nodes], :xyz)...)'




        #Lower triangle on super-element defined by overlapMat()
        for jj = 1:size(R_array, 3)

            #use each kernel region as a clipping polygon. The reference polygon is the element
            elmtInt = sutherlandHodgman(tL[:, 1:2], R_array[:, :, jj])
            #If vertices are near equivalent, remove duplicates
            if ~isempty(elmtInt) && size(elmtInt, 1) > 3
                elmtInt = reducePolyJ(elmtInt, 10^(-14))
            end
            #Only perform integrations if the integral region is two-dimensional
            if size(elmtInt, 1) > 2

                if ~isempty(elmtInt)
                    #If the quadrature region has more than 3 vertices, break it up into triangles and integrate over those
                    if size(elmtInt, 1) > 3
                        Tess = Fan_triangulation(elmtInt)
                        #Tess=Del_triangulation(elmtInt);
                        for z = 1:size(Tess, 3)

                            #Map to global domain
                            XX, YY = cartesian_to_global(xi1, xi2, Tess[:, :, z])
                            Jac = Jacobian_Tri(Tess[:, :, z])

                            #Map quadrature nodes to reference triangle of the element
                            XI1, XI2 = global_to_cartesian(XX, YY, tL)
                            # Now map to the reference square
                            ETA1, ETA2 = cartesian_to_collapsed(XI1, XI2)


							uhatVal = phys_at_element(fld,ETA1, ETA2,KL)
                            #uhatVal = evaluate_Dubiner_modal(p, data.modes[:, KL], ETA1, ETA2)
                            kernelVal = evalKernelPeriodic(kernels, Hx, Hy, XX, YY, offset, Lx, Ly, loc_array[jj], R_coords)
							#@info "kernel $kernelVal, weights $wl, $wr, Jacobian $Jac"
                            intVal = intVal + Jac * sum(wl[zx] .* wr[zy] .* uhatVal .* kernelVal .* (0.5 * (1 .- zr[zy])))
                        end
                    else
                        #Map to global domain
                        XX, YY = cartesian_to_global(xi1, xi2, elmtInt[:, :])
                        Jac = Jacobian_Tri(elmtInt[:, :])
                        #Map quadrature nodes to reference triangle of the element
                        XI1, XI2 = global_to_cartesian(XX, YY, tL)
                        # Now map to the reference square
                        ETA1, ETA2 = cartesian_to_collapsed(XI1, XI2)

                        uhatVal = phys_at_element(fld,ETA1, ETA2,KL)
                        #uhatVal = evaluate_Dubiner_modal(p, data.modes[:, KL], ETA1, ETA2)

                        kernelVal = evalKernelPeriodic(kernels, Hx, Hy, XX, YY, offset, Lx, Ly, loc_array[jj], R_coords)
                        intVal = intVal + Jac * sum(wl[zx] .* wr[zy] .* uhatVal .* kernelVal .* (0.5 * (1 .- zr[zy])))


                    end
                end
            end
        end

        #upper triangle on super-element defined by overlapMat()
        for jj = 1:size(R_array, 3)

            elmtInt = sutherlandHodgman(tR[:, 1:2], R_array[:, :, jj])
            if ~isempty(elmtInt) && size(elmtInt, 1) > 3
                elmtInt = reducePolyJ(elmtInt, 10^(-14))
            end
            if size(elmtInt, 1) > 2
                if ~isempty(elmtInt)

                    if size(elmtInt, 1) > 3
                        Tess = Fan_triangulation(elmtInt)
                        #Tess=Del_triangulation(elmtInt);
                        for z = 1:size(Tess, 3)
                            #Map to global domain
                            XX, YY = cartesian_to_global(xi1, xi2, Tess[:, :, z])
                            Jac = Jacobian_Tri(Tess[:, :, z])

                            #Map quadrature nodes to reference triangle of the element
                            XI1, XI2 = global_to_cartesian(XX, YY, tR)
                            # Now map to the reference square
                            ETA1, ETA2 = cartesian_to_collapsed(XI1, XI2)


							uhatVal = phys_at_element(fld,ETA1, ETA2,KR)
                            #uhatVal = evaluate_Dubiner_modal(p, data.modes[:, KR], ETA1, ETA2)

                            kernelVal = evalKernelPeriodic(kernels, Hx, Hy, XX, YY, offset, Lx, Ly, loc_array[jj], R_coords)
                           intVal = intVal + Jac * sum(wl[zx] .* wr[zy] .* uhatVal .* kernelVal .* (0.5 * (1 .- zr[zy])))



                        end
                    else

                        #Map to global domain

                        XX, YY = cartesian_to_global(xi1, xi2, elmtInt[:, :])
                        Jac = Jacobian_Tri(elmtInt[:, :])

                        #Map quadrature nodes to reference triangle of the element
                        XI1, XI2 = global_to_cartesian(XX, YY, tR)
                        # Now map to the reference square
                        ETA1, ETA2 = cartesian_to_collapsed(XI1, XI2)


                        uhatVal = phys_at_element(fld,ETA1, ETA2,KR)
                        #uhatVal = evaluate_Dubiner_modal(p, data.modes[:, KR], ETA1, ETA2)

                        kernelVal = evalKernelPeriodic(kernels, Hx, Hy, XX, YY, offset, Lx, Ly, loc_array[jj], R_coords)
                        intVal = intVal + Jac * sum(wl[zx] .* wr[zy] .* uhatVal .* kernelVal .* (0.5 * (1 .- zr[zy])))


                    end
                end
            end

        end
    end


    ##################################################


    # rescale the integral by the kernel scaling
    out = intVal / (Hx * Hy)

    return out
end


#Temporary
function cartesian_to_collapsed(xi1,xi2)
    #Maps triangular reference element to
    # collapsed coordinate system's reference quadrilateral

    eta1=2 .* (1 .+ xi1) ./ (1 .- xi2) .- 1
    eta2=xi2

    return eta1,eta2
end



function cartesian_to_global(xi1,xi2,vert)

    #Map reference quadrilateral to global coordinates
    vx1=vert[1,1]
    vx2=vert[2,1]
    vx3=vert[3,1]
    vy1=vert[1,2]
    vy2=vert[2,2]
    vy3=vert[3,2]



    x=0.5 .* (-(xi1 .+ xi2) .* vx1 .+ (xi1 .+ 1) .* vx2 .+ (xi2 .+ 1) .* vx3)
    y=0.5 .* (-(xi1 .+ xi2) .* vy1 .+ (xi1 .+ 1) .* vy2 .+ (xi2 .+ 1) .* vy3)
    return x,y
end

function collapsed_to_cartesian(eta1,eta2)
    #Maps collapsed coordinate system's reference quadrilateral to the
    #triangular reference element

    xi2=eta2
    xi1=0.5 .* (1 .+ eta1) .* (1 .- eta2) .- 1

    return xi1,xi2
end


function global_to_cartesian(x,y,vert)

    #Map reference quadrilateral to global coordinates
    vx1=vert[1,1]
    vx2=vert[2,1]
    vx3=vert[3,1]
    vy1=vert[1,2]
    vy2=vert[2,2]
    vy3=vert[3,2]

    det=(vx2 .- vx1) .* (vy3 .- vy1) .- (vx3 .- vx1) .* (vy2 .- vy1)


    xi1=1 ./ det .* ((2 .* x .- vx2 .- vx3) .* (vy3 .- vy1) .+ (2 .* y .- vy2 .- vy3) .* (vx1 .- vx3))
    xi2=1 ./ det .* ((2 .* x .- vx2 .- vx3) .* (vy1 .- vy2) .+ (2 .* y .- vy2 .- vy3) .* (vx2 .- vx1))

    return xi1,xi2
end
function funcBreakPoints_periodic(boundaries, support)
    # Returns a list (in world coords) of the periodic breakpoints of the
    # function as defined by the boundaries and the region of support
    # doesn't include the actual support limits -- strictly internal breaks
    #
    # assumes that the size of the support is less than that of the domain

    breaks = []
    lower = support[1]
    upper = support[2]


    domainMin = boundaries[1]
    domainMax = boundaries[length(boundaries)]
    domainWidth = domainMax - domainMin

    indicator = 1
    while domainMax * indicator < upper
        indicator = indicator + 1
    end
    indicator = 2 * indicator
    #     We need our extension to contain the support. The indicator states
    #     how many time we need to add the original domainWidth to contain the
    #     upper end of the support. As I assumed the support is symmetric about
    #     the domain, doubling this will enable our total extended interval
    #     width to contain the whole support. We want to extend the tesselation
    #     of the domain out through the upper bound of the support. To do so we
    #     periodically extend our boundary nodes as follows
    Boundmatrix = zeros(indicator, length(boundaries) - 1)
    for i = 1:indicator
        Boundmatrix[i, :] = boundaries[2:length(boundaries)] .+ i .* domainWidth
    end
    Boundmatrix2 = Boundmatrix'
    extension = vec(Boundmatrix2[:]')
    doubleBoundary = [boundaries; extension]
    # shift support's lower end into the domain
    offset = 0     # the number of support widths I had to shift the support to get it to start in the domain
    while (lower < domainMin)
        offset = offset + 1
        lower = lower + domainWidth
        upper = upper + domainWidth
    end
    while (lower > domainMax)
        offset = offset - 1
        lower = lower - domainWidth
        upper = upper - domainWidth
    end



    boundaryIndex = 1
    while (doubleBoundary[boundaryIndex] <= lower)
        boundaryIndex = boundaryIndex + 1
    end


    brkIndex = 1
    while (doubleBoundary[boundaryIndex] < upper)
        breaks = append!(breaks, doubleBoundary[boundaryIndex])
        boundaryIndex = boundaryIndex + 1
        brkIndex = brkIndex + 1
    end



    breaks = breaks .- domainWidth .* offset

    return breaks
    #    breaks = [breaks(1)-supportWidth/(size(breaks,2)) breaks breaks(size(breaks,2))+supportWidth/(size(breaks,2))];
end

function periodic_mesh_shift(boundaries, brks)


    domainMin = minimum(boundaries)
    domainMax = maximum(boundaries)
    domainWidth = domainMax - domainMin
    out = zeros(size(brks))

    for i = 1:length(brks)
        out[i] = brks[i]
        # map into computational domain
        while (out(i) < domainMin)
            out[i] = out[i] + domainWidth
        end
        while (out[i] > domainMax)
            out[i] = out[i] - domainWidth
        end
    end

    return out
end


function Jacobian_Tri(vert)

    x_xi1 = 0.5 * (vert[2, 1] - vert[1, 1])
    x_xi2 = 0.5 * (vert[3, 1] - vert[1, 1])
    y_xi1 = 0.5 * (vert[2, 2] - vert[1, 2])
    y_xi2 = 0.5 * (vert[3, 2] - vert[1, 2])

    J = -x_xi2 .* y_xi1 + x_xi1 .* y_xi2

    return J
end



function kspptPeriodic(xmin, xmax, ymin, ymax, VNW, VSW, VSE, VNE)

    #Function that identifies the no longer connected regions on the global
    #mesh under the kernel support for periodic boundary conditions.
    # Note that all this is doing is applying the clipping algorithms to
    # translates of the kernel support.

    Lx = xmax - xmin
    Ly = ymax - ymin

    shift_array = [0 0
        1 0
        1 1
        0 1
        -1 0
        0 -1
        -1 -1
        -1 1
        1 -1]

    domain_poly = [xmin ymin
        xmax ymin
        xmax ymax
        xmin ymax]

    Kernel_poly = [VSW; VSE; VNE; VNW]

    R_array = []
    for ind = 1:9
        #Need to re-sort after shifting
        clip_poly = putCCorder(Kernel_poly .+ ([Lx; Ly] .* shift_array[ind, :])' .* ones(size(Kernel_poly)))


        R_proxy = sutherlandHodgman(domain_poly, clip_poly)

        if isempty(R_proxy)
        elseif isempty(R_array)
            R_array = R_proxy
        else
            R_array = cat([R_array, R_proxy]..., dims=3)
        end
    end


    return R_array
end



function superElmtToTri(I, J, Nx)
    # Function that when given super-element ordering computes the single index
    # for the two triangular sub-elements  composing that quadrilateral
    # super-element.

    KL = (2 * I - 1) + (J - 1) * (2 * Nx)
    KR = (2 * I) + (J - 1) * (2 * Nx)

    return KL, KR
end


function kFootPrint2(kHalfWidth_x, kHalfWidth_y, boundaries_x, boundaries_y, offset_x, offset_y)
    #Routine return the (I,J) indices of the quadrilateral super-element
    #overlapped by the kernel support. Note this assumes a uniform mesh, but
    #not a uniform kernel scaling in the x and y directions.
    #Note that owing to the mesh uniformity, I can just compute the left and
    #bottom most element and the top right most element and just vary the
    #indices
    #Breaks if kernel support larger than the domain (L<(RS+BS-1)H)

    #Put kernel sppt into computational domain

    Nx = length(boundaries_x) - 1
    Ny = length(boundaries_y) - 1

    x_left = offset_x - kHalfWidth_x
    y_bottom = offset_y - kHalfWidth_y
    x_right = offset_x + kHalfWidth_x
    y_top = offset_y + kHalfWidth_y

    Lx = boundaries_x[end] - boundaries_x[1]
    Ly = boundaries_y[end] - boundaries_y[1]

    I_left, J_bottom = findQuadSuperElmntPeriodic(x_left, y_bottom, boundaries_x, boundaries_y)
    I_right, J_top = findQuadSuperElmntPeriodic(x_right, y_top, boundaries_x, boundaries_y)


    I_vec = [I_left:1:I_right;]
    J_vec = [J_bottom:1:J_top;]

    #Meshgrid wont work as I want to store the individual element indices in a
    #M x 2 matrix
    count = 1
    Index_mat = zeros(length(I_vec) * length(J_vec), 2)
    for i = 1:length(I_vec)

        if I_vec[i] > Nx
            I_vec[i] = I_vec[i] - Nx
        elseif I_vec[i] < 1
            I_vec[i] = I_vec[i] + Nx
        end

        for j = 1:length(J_vec)

            if J_vec[j] > Ny
                J_vec[j] = J_vec[j] - Ny
            elseif J_vec[j] < 1
                J_vec[j] = J_vec[j] + Ny
            end



            Index_mat[count, 1] = I_vec[i]
            Index_mat[count, 2] = J_vec[j]
            count = count + 1
        end
    end

    return Index_mat

end


function findQuadSuperElmntPeriodic(x, y, boundaries_x, boundaries_y)
    #Determine what quadrilateral super-element (x,y) belongs to. Supposes the
    #domain need only be extended once in either either direction

    Nx = length(boundaries_x) - 1
    Ny = length(boundaries_y) - 1
    Lx = boundaries_x[end] - boundaries_x[1]
    Ly = boundaries_y[end] - boundaries_y[1]

    bndsx = [boundaries_x[1:end-1] .- Lx; boundaries_x[1:end-1]; boundaries_x .+ Lx]
    bndsy = [boundaries_y[1:end-1] .- Ly; boundaries_y[1:end-1]; boundaries_y .+ Ly]

    I = 1 - Nx
    J = 1 - Ny
    cntx = 1
    cnty = 1
    while x > bndsx[cntx+1]
        I = I + 1
        cntx = cntx + 1
    end

    while y > bndsy[cnty+1]
        J = J + 1
        cnty = cnty + 1
    end

    return I, J
end



function kFootPrint(I, J, kHalfWidth_x, kHalfWidth_y, hx, hy, Nx, Ny)
    #Routine return the (I,J) indices of the quadrilateral super-element
    #overlapped by the kernel support. Note this assumes a uniform mesh, but
    #not a uniform kernel scaling in the x and y directions.
    #Note that owing to the mesh uniformity, I can just compute the left and
    #bottom most element and the top right most element and just vary the
    #indices

    I_left = (I) - ceil(kHalfWidth_x / hx)
    I_right = I + ceil(kHalfWidth_x / hx)
    J_bottom = (J) - ceil(kHalfWidth_y / hy)
    J_top = J + ceil(kHalfWidth_y / hy)

    I_vec = [I_left:1:I_right;]
    J_vec = [J_bottom:1:J_top;]
    #Meshgrid wont work as I want to store the individual element indices in a
    #M x 2 matrix
    Index_mat = zeros(length(I_vec) * length(J_vec), 2)

    count = 1
    for i = 1:length(I_vec)
        for j = 1:length(J_vec)
            Index_mat[count, 1] = I_vec[i]
            Index_mat[count, 2] = J_vec[j]
            count = count + 1
        end
    end




    A = size(Index_mat, 1)
    for a = 1:A

        x_ind = Index_mat[a, 1]
        y_ind = Index_mat[a, 2]
        # map into computational domain
        while (x_ind < 1)
            x_ind = x_ind + Nx
        end
        while (x_ind > Nx)
            x_ind = x_ind - Nx
        end

        while (y_ind < 1)
            y_ind = y_ind + Ny
        end
        while (y_ind > Ny)
            y_ind = y_ind - Ny
        end

        Index_mat[a, 1] = x_ind
        Index_mat[a, 2] = y_ind
    end


    return Index_mat

end


function putCCorder(poly)
    x = poly[:, 1]
    y = poly[:, 2]
    cx = mean(x)
    cy = mean(y)
    a = atan.(y .- cy, x .- cx)
    order = sortperm(a)
    x = x[order]
    y = y[order]
    out = zeros(length(x), 2)
    out[:, 1] = x
    out[:, 2] = y

    return out
end

function reducePolyJ(Poly, tol)
    #Given polynomial vertices and a tolerance parameter, removes vertices that
    #are duplicates up to that tolerance. The MATLAB version was better and used dynamic programming :(


    N = size(Poly, 1)
    X_list = Poly[1, 1]
    Y_list = Poly[1, 2]
    for i = 2:N
        p_x = []
        p_y = []
        for j = 1:size(X_list, 1)
            p_x = []
            p_y = []
            if (abs(X_list[j] - Poly[i, 1]) <= tol) && (abs(Y_list[j] - Poly[i, 2]) <= tol)
                break
            end
            p_x = Poly[i, 1]
            p_y = Poly[i, 2]

        end
        if ~isempty(p_x)
            X_list = hcat([X_list, p_x]...)
            Y_list = hcat([Y_list, p_y]...)
        end
    end
    verifiedVert = [X_list' Y_list']

    return verifiedVert
end


function findQuadSuperElmnt(x, y, boundaries_x, boundaries_y)

    #Determine what quadrilateral super-element (x,y) belongs to. Assumes
    #semi-uniform mesh and wraps periodically.

    Nx = length(boundaries_x) - 1
    Ny = length(boundaries_y) - 1

    Lx = boundaries_x[end] - boundaries_x[1]
    Ly = boundaries_y[end] - boundaries_y[1]

    I = 1
    while x > boundaries_x[end]
        x = x - Lx
        I = I + Nx
    end
    while x < boundaries_x[1]
        x = x + Lx
        I = I - Nx
    end
    cx = 1
    while x > boundaries_x[cx+1]
        I = I + 1
        cx = cx + 1
    end

    J = 1
    while y > boundaries_y[end]
        y = y - Ly
        J = J + Ny
    end
    while y < boundaries_y[1]
        y = y + Ly
        J = J - Ny
    end
    cy = 1
    while y > boundaries_y[cy+1]
        J = J + 1
        cy = cy + 1
    end


    return I, J

end


function Fan_triangulation(V)
    #Function that takes N*2 matrix of vertices on an N-gon and returns a
    #tensor array containing the vertices of N-2 triangles tesselating that
    #N-gon. This tesselation is not unique and it is only guaranteed to work for convex
    #polygons

    N = size(V, 1)

    Tess = zeros(3, 2, N - 2)
    for i = 1:N-2
        Tess[1, :, i] = V[1, :]
        Tess[2, :, i] = V[i+1, :]
        Tess[3, :, i] = V[i+2, :]
    end

    return Tess
end


function evalKernelPeriodic(kernels, Hx, Hy, x, y, offset, Lx, Ly, R_loc, R_coords)
    #Evaluates the kernel periodically, but supposes all the eval points either
    #need to be or do not need to be periodically shifted. This is ok because
    #it should only be called over a single element.

    tol = 10^(-14)
    #tol=0

    #Makes the output an array of the same size as x
    out = zeros(size(x))

    VNW = R_coords[1, :, R_loc]
    VSW = R_coords[2, :, R_loc]
    VSE = R_coords[3, :, R_loc]
    VNE = R_coords[4, :, R_loc]

    while minimum(x) < VNW[1] - tol
        x = x .+ Lx
    end

    while maximum(x) > VNE[1] + tol
        x = x .- Lx
    end

    while minimum(y) < VSW[2] - tol
        y = y .+ Ly
    end

    while maximum(y) > VNW[2] + tol
        y = y .- Ly
    end

    xinput=round.(((x .- offset[1]) ./ Hx), digits=14)
    yinput=round.(((y .- offset[2]) ./ Hy), digits=14)
    for i=1:length(x)
        if xinput[i]==-0.0
            xinput[i]=0.0
        end
        if yinput[i]==-0.0
            yinput[i]=0.0
        end
    end
    out = evaluate_kernel(kernels[1].T, kernels[1].c, xinput) .* evaluate_kernel(kernels[2].T, kernels[2].c, yinput)
    return out
end



#The inputs are a table of x-y pairs for the verticies of the subject
#polygon and boundary polygon. (x values in column 1 and y values in column
#2) The output is a table of x-y pairs for the clipped version of the
#subject polygon.

## This is based on a MATLAB algorithm found at http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#MATLAB_.2F_Octave,
#though I have editted it to use more robust intersection algorithms, and to deal with Julia vector/matrix conventions.



function sutherlandHodgman(subjectPolygon, clipPolygon)

    ## Sutherland-Hodgman Algorithm
    clipPolygon_x = clipPolygon[:, 1]
    clipPolygon_y = clipPolygon[:, 2]

    clippedPolygon_x = subjectPolygon[:, 1]
    clippedPolygon_y = subjectPolygon[:, 2]

    numVerticies = size(clipPolygon, 1)
    clipVertexPrevious_x = clipPolygon[end, 1]
    clipVertexPrevious_y = clipPolygon[end, 2]

    for clipVertex = 1:numVerticies

        clipBoundary_x = vcat([clipPolygon_x[clipVertex], clipVertexPrevious_x]...)
        clipBoundary_y = vcat([clipPolygon_y[clipVertex], clipVertexPrevious_y]...)

        inputList_x = clippedPolygon_x
        inputList_y = clippedPolygon_y

        clippedPolygon_x = []
        clippedPolygon_y = []
        if ~isempty(inputList_x)
            previousVertex_x = inputList_x[end]
            previousVertex_y = inputList_y[end]
        end

        for subjectVertex = 1:length(inputList_x)

            if inside(inputList_x[subjectVertex], inputList_y[subjectVertex], clipBoundary_x, clipBoundary_y) == 1

                if ~(inside(previousVertex_x, previousVertex_y, clipBoundary_x, clipBoundary_y)) == 1
                    subjectLineSegment_x = vcat([previousVertex_x, inputList_x[subjectVertex]]...)
                    subjectLineSegment_y = vcat([previousVertex_y, inputList_y[subjectVertex]]...)
                    subjectLineSegment = zeros(2, 2)
                    subjectLineSegment[:, 1] = subjectLineSegment_x
                    subjectLineSegment[:, 2] = subjectLineSegment_y
                    clipBoundary = hcat([clipBoundary_x, clipBoundary_y]...)

                    edge_int = LineSegmentIntersection(clipBoundary, subjectLineSegment)


                    clippedPolygon_x = cat([clippedPolygon_x, edge_int[1]]..., dims=1)
                    clippedPolygon_y = cat([clippedPolygon_y, edge_int[2]]..., dims=1)

                end


                clippedPolygon_x = vcat([clippedPolygon_x, inputList_x[subjectVertex]]...)
                clippedPolygon_y = vcat([clippedPolygon_y, inputList_y[subjectVertex]]...)


            elseif inside(previousVertex_x, previousVertex_y, clipBoundary_x, clipBoundary_y) == 1
                subjectLineSegment_x = vcat([previousVertex_x, inputList_x[subjectVertex]]...)
                subjectLineSegment_y = vcat([previousVertex_y, inputList_y[subjectVertex]]...)
                subjectLineSegment = zeros(2, 2)

                subjectLineSegment[:, 1] = subjectLineSegment_x
                subjectLineSegment[:, 2] = subjectLineSegment_y
                clipBoundary = hcat([clipBoundary_x, clipBoundary_y]...)

                edge_int = LineSegmentIntersection(clipBoundary, subjectLineSegment)


                clippedPolygon_x = cat([clippedPolygon_x, edge_int[1]]..., dims=1)
                clippedPolygon_y = cat([clippedPolygon_y, edge_int[2]]..., dims=1)


            end

            previousVertex_x = inputList_x[subjectVertex]
            previousVertex_y = inputList_y[subjectVertex]

            clipVertexPrevious_x = clipPolygon_x[clipVertex]
            clipVertexPrevious_y = clipPolygon_y[clipVertex]

        end #for subject verticies
    end #for boundary verticies

    clippedPolygon = hcat([clippedPolygon_x, clippedPolygon_y]...)
    return clippedPolygon
end #sutherlandHodgman


function inside(point_x, point_y, boundary_x, boundary_y)

    indicator = (boundary_x[2] - boundary_x[1]) * (point_y - boundary_y[1]) - (boundary_y[2] - boundary_y[1]) * (point_x - boundary_x[1])

    if (indicator <= 0)
        in = true
    else
        in = false
    end
    return in

end

function LineSegmentIntersection(L1, L2)
    #2D only. Supposes that the lines defined by the line-segements are not
    #collinear [In which case you only get 1 intersection] or parallel [In
    #which case there is no solution]

    intersection = zeros(1, 2)


    #intersection via parameterization


    t = (((L1[1, 1] - L2[1, 1]) * (L2[1, 2] - L2[2, 2]) - (L1[1, 2] - L2[1, 2]) * (L2[1, 1] - L2[2, 1])) /
         ((L1[1, 1] - L1[2, 1]) * (L2[1, 2] - L2[2, 2]) - (L1[1, 2] - L1[2, 2]) * (L2[1, 1] - L2[2, 1])))

    u = (((L1[1, 1] - L2[1, 1]) * (L1[1, 2] - L1[2, 2]) - (L1[1, 2] - L2[1, 2]) * (L1[1, 1] - L1[2, 1])) /
         ((L1[1, 1] - L1[2, 1]) * (L2[1, 2] - L2[2, 2]) - (L1[1, 2] - L1[2, 2]) * (L2[1, 1] - L2[2, 1])))

    if (u > 1) || (u < 0) || (t > 1) || (t < 0)
        #     fprintf['Segments do not intersect \n']
        #     pause
        # return []
    end
    intersection[1] = L1[1, 1] + t * (L1[2, 1] - L1[1, 1])
    intersection[2] = L1[1, 2] + t * (L1[2, 2] - L1[1, 2])

    #Equivalent answer depending on the second parameterized line
    # intersection2[1] = L2[1,1]+u*(L2[2,1]-L2[1,1])
    # intersection2[2] = L2[1,2]+u*(L2[2,2]-L2[1,2])

    return intersection
end



function standard_coords_stt(z, v)

    d, n = size(v)
    mz = 0.5 .* (1.0 .- z)
    pz = 0.5 .* (1.0 .+ z)
    if n == 4 && d == 3
        return v[:, 1] .* mz[1] .* mz[2] .+ v[:, 2] .* pz[1] .* mz[2] .+
               v[:, 4] .* mz[1] .* pz[2] .+ v[:, 3] .* pz[1] .* pz[2]
    end
    if n == 3 && d == 3

        #First map collapsed to cartesian
        xi2 = z[2]
        xi1 = 0.5 * (1 + z[1]) .* (1 - z[2]) - 1
        #Map reference quadrilateral to global coordinates
        vx1 = v[1, 1]
        vx2 = v[1, 2]
        vx3 = v[1, 3]
        vy1 = v[2, 1]
        vy2 = v[2, 2]
        vy3 = v[2, 3]

        x = 0.5 * (-(xi1 + xi2) .* vx1 + (xi1 + 1) .* vx2 + (xi2 + 1) .* vx3)
        y = 0.5 * (-(xi1 + xi2) .* vy1 + (xi1 + 1) .* vy2 + (xi2 + 1) .* vy3)
        return [x, y]
    end
end
