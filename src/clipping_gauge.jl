"""
    check_preclip(t,n)
For an n-qubit stabilizer state specified by the tableau,
check that the preclip condition is satisfied.
"""
function check_preclip(my_tableau, n)
    stab = my_tableau[(n+1):(2*n),1:(2n)]
    return check_preclip_stab(stab, n)
end

"""
    check_preclip_stab(stab,n)
For an n-qubit stabilizer state specified by the stabilizer
part of the tableau, check that the preclip condition is satisfied.
"""
function check_preclip_stab(stab, n)
    num_stab = size(stab)[1]
    xAll = stab[:, 1:n]
    zAll = stab[:, (n+1):(2*n)]
    vals = xAll .| zAll
    lefts = [findfirst(vals[i,:]) for i = 1:num_stab]
    for i = 1:n
        c = count(lefts .== i)
        if c > 2
            return false
        elseif c == 2
            inds = findall(lefts .== i)
            if (xAll[inds[1],i] == xAll[inds[2],i]) && (zAll[inds[1],i] == zAll[inds[2],i])
                return false
            end
        end
    end
    return true
end

"""
    check_preclip_stab_right(stab,n)
For an n-qubit stabilizer state specified by the stabilizer part
of the tableau, check that the right-hand preclip condition is satisfied.
"""
function check_preclip_stab_right(stab, n)
    num_stabs = size(stab)[1]
    xAll = stab[:, 1:n]
    zAll = stab[:, (n+1):(2*n)]
    vals = xAll .| zAll
    rights = [findlast(vals[i,:]) for i = 1:num_stabs]
    for i = 1:n
        c = count(rights .== i)
        if c > 2
            return false
        elseif c == 2
            inds = findall(rights .== i)
            if (xAll[inds[1],i] == xAll[inds[2],i]) && (zAll[inds[1],i] == zAll[inds[2],i])
                return false
            end
        end
    end
    return true
end

"""
    check_clip(t,n)
For an n-qubit stabilizer state specified by the tableau,
check that the clipped gauge condition is satisfied.
"""
function check_clip(my_tableau, n; verbose=false)
    stab = my_tableau[(n+1):(2*n), 1:(2*n)]
    return check_clip_stab(stab, n, verbose=verbose)
end

"""
    check_clip_stab(stab,n)
For an n-qubit stabilizer state specified by the stabilizer part
of the tableau, check that the clipped gauge condition is satisfied.
"""
function check_clip_stab(stab, n; verbose=false)
    num_stabs = size(stab)[1]
    flag = true
    xAll = stab[:,1:n]
    zAll = stab[:,(n+1):(2*n)]
    vals = xAll .| zAll
    lefts = [findfirst(vals[i,:]) for i = 1:num_stabs]
    rights = [findlast(vals[i,:]) for i = 1:num_stabs]
    for i = 1:n
        cL = count(lefts .== i)
        cR = count(rights .== i)
        if (cL + cR) != 2
            if verbose
                println("cL + cR = ", cL+cR, " at index ", i)
            end
            flag = false
            # return false
        end
        if cL == 2
            inds = findall(lefts .== i)
            if (xAll[inds[1],i] == xAll[inds[2],i]) && (zAll[inds[1],i] == zAll[inds[2],i])
                if verbose
                    println("both cL same at index ", i)
                end
                flag = false
                # return false
            end
        end
        if cR == 2
            inds = findall(rights .== i)
            if (xAll[inds[1],i] == xAll[inds[2],i]) && (zAll[inds[1],i] == zAll[inds[2],i])
                if verbose
                    println("both cR same at index ", i)
                end
                flag = false
                # return false
            end
        end
    end
    return flag
    # return true
end

"""
    Z2_gaussian_eliminate(stab,n)
For an n-qubit stabilizer state specified by the stabilizer
part of the tableau, perform Gaussian elimination, proceeding
from left to right. Afterward, stab should satisfy the left
preclip condition.
"""
function Z2_gaussian_eliminate(stab,n)
    num_stabs = size(stab)[1]
    M = falses(num_stabs,2*n)
    M[:,1:2:2*n] = stab[:,1:n]
    M[:,2:2:2*n] = stab[:,(n+1):(2*n)]

    i = 1
    for j = 1:2*n
        p = findfirst(M[i:num_stabs, j])
        if p != nothing
            p = p + i - 1
            tmp = copy(M[i,:])
            M[i, :] .= M[p, :]
            M[p, :] .= tmp
            for k = 1:num_stabs
                if M[k,j] == 1 && i != k
                    M[k,:] .= (M[k,:] .⊻ M[i,:])
                end
            end
            i += 1
        end
    end

    new_M = falses(num_stabs, 2*n)
    for i=1:num_stabs
        new_op = falses(2*n)
        new_op[1:n] .= M[i,1:2:2*n]
        new_op[(n+1):(2*n)] .= M[i,2:2:2*n]
        new_M[i,:] .= new_op
    end
    @assert check_preclip_stab(new_M, n)
    return new_M
end

"""
    clipping_gauge(stab,n)
For an n-qubit stabilizer state specified by the stabilizer
part of the tableau, go to the clipping gauge by performing
Gaussian elimination from right to left, always eliminating
the longer stabilizer by the shorter stabilizer.
"""
function clipping_gauge(stab, n)
    num_stabs = size(stab)[1]
    @assert num_stabs <= n
    if num_stabs <= 1
        return stab
    end

    new_stab = Z2_gaussian_eliminate(stab, n)
    @assert num_stabs == size(new_stab)[1]

    M = falses(num_stabs,2*n)
    M[:,1:2:2*n] = new_stab[:,1:n]
    M[:,2:2:2*n] = new_stab[:,(n+1):(2*n)]

    i = 1
    for j = 2n:-1:1
        if i > num_stabs
            break
        end

        s = sum(M,dims=2)
        if any(s .== 0)
            throw(ErrorException("At least one row is empty"))
        end

        p = findall(M[i:num_stabs,j])
        if length(p) == 0
            continue
        end
        if length(p) == 1
            p = p[1] + i - 1
            tmp = copy(M[i,:])
            M[i,:] .= M[p, :]
            M[p, :] .= tmp
            i += 1
            continue
        end

        p = p .+ (i-1)

        # find shortest
        indShortest = -1
        lMax = -1
        for ind in p
            l = findfirst(M[ind,:])
            if l > lMax
                indShortest = ind
                lMax = l
            end
        end
        if indShortest < 1
            throw(ErrorException("Invalid index for shortest stabilizer"))
        end

        tmp = copy(M[i,:])
        M[i,:] .= M[indShortest, :]
        M[indShortest, :] .= tmp

        for k = (i+1):num_stabs
            if M[k,j] == 1 && i != k
                if M[k,:] == M[i,:]
                    throw(ErrorException("Eliminating repeated row"))
                end
                M[k,:] .= (M[k,:] .⊻ M[i,:])
            end
        end

        i += 1

    end
    new_stab = falses(num_stabs, 2*n)
    for i=1:num_stabs
        new_op = falses(2*n)
        new_op[1:n] .= M[i,1:2:2*n]
        new_op[(n+1):(2*n)] .= M[i,2:2:2*n]
        new_stab[i,:] .= new_op
    end
    @assert check_preclip_stab(new_stab, n)
    @assert check_preclip_stab_right(new_stab, n)
    @assert check_clip_stab(new_stab, n) | (num_stabs < n)
    return new_stab
end

"""
    tableau_to_stab_lengths(t,n)
For an n-qubit stabilizer state specified by the tableau,
convert to the clipping gauge and then compute the stabilizer
lengths.

Optional parameter pure (boolean) determines whether we are returning
the lengths of the stabilizers (pure) or the generators of the mixed space.
"""
function tableau_to_stab_lengths(t, n; pure=true)

    stab = copy(t[(n+1):(2n), :])
    if pure
        mask = (stab[:,2n+2] .== 0)
    else
        mask = (stab[:,2n+2] .== 1)
    end
    stab = stab[mask, 1:(2n)]
    num_stabs = size(stab)[1]

    stab = clipping_gauge(stab, n)

    xAll = stab[:, 1:n]
    zAll = stab[:, (n+1):(2*n)]
    vals = xAll .| zAll
    lefts = [findfirst(vals[i,:]) for i = 1:num_stabs]
    rights = [findlast(vals[i,:]) for i = 1:num_stabs]

    stab_lens = rights .- lefts .+ 1
    return stab_lens
end

"""
    tableau_to_stab_ends(t,n)
For an n-qubit stabilizer state specified by the tableau, convert to the
clipping gauge and then compute the density of stabilizer endpoints.
"""
function tableau_to_stab_ends(t, n)

    stab = copy(t[(n+1):(2n), :])
    is_pure = (stab[:,2n+2] .== 0)
    stab = stab[is_pure, 1:(2n)]
    num_stabs = size(stab)[1]

    stab = clipping_gauge(stab, n)

    xAll = stab[:, 1:n]
    zAll = stab[:, (n+1):(2*n)]
    vals = xAll .| zAll
    lefts = [findfirst(vals[i,:]) for i = 1:num_stabs]
    rights = [findlast(vals[i,:]) for i = 1:num_stabs]

    endpoints_density = [count(lefts .== i) + count(rights .== i) for i = 1:n]
    return endpoint_density
end
