using Random, LinearAlgebra, Statistics
using DelimitedFiles, ArgParse, Printf
using Plots
ENV["GKSwstype"]="nul"
using LoopVectorization

function test_func_2(x)
    return x + 3
end

function g(x1,z1,x2,z2)
    # Helper function for rowsum
    if x1 == z1 == 0
        return 0
    elseif x1 == z1 == 1
        return z2 - x2
    elseif x1 == 1 && z1 == 0
        return z2 * (2 * x2 - 1)
    else
        return x2 * (1 - 2 * z2)
    end
end

@inline function rowsum(tableau::BitArray, h::Int, i::Int, n::Int)
    # Set rh the phase register
    gSum = mod(sum(g.(tableau[i,1:n],tableau[i,(n+1):2n],tableau[h,1:n],tableau[h,(n+1):2n])), 4)
    tableau[h,2n+1] = (mod(gSum + 2*(tableau[h,2n+1] + tableau[i,2n+1]), 4) != 0)
    # Set row entries to xor
    tableau[h,1:(2n)] .= (tableau[i,1:(2n)] .⊻ tableau[h,1:(2n)])
    return tableau
end

function check_commute(op1::BitArray, op2::BitArray, n::Int)

    gSum1 = mod(sum(g.(op1[1:n],op1[(n+1):2n],op2[1:n],op2[(n+1):2n])),4)
    gSum2 = mod(sum(g.(op2[1:n],op2[(n+1):2n],op1[1:n],op1[(n+1):2n])),4)

    return gSum1 == gSum2

end

function print_operator(stab, n)
    str = ""
    for j = 1:n
        if stab[j] == 1 && stab[j+n] == 1
            str = string(str, "Y")
        elseif stab[j] == 1 && stab[j+n] == 0
            str = string(str, "X")
        elseif stab[j+n] == 1
            str = string(str,"Z")
        else
            str = string(str,"I")
        end
    end
    println(str)
end

function measure_operator(tableau::BitArray, op::BitArray, n::Int)
    p = 0
    @inbounds for i = (n+1):(2n)
        if check_commute(tableau[i,:], op, n) == false
            p = i
            break
        end
    end

    if p > 0
        @inbounds for i = 1:2n
            if (i != p) && (check_commute(tableau[i,:], op, n) == false)
                tableau = rowsum(tableau, i, p, n)
            end
        end
        tableau[p-n, :] .= copy(tableau[p,:])
        tableau[p, :] .= copy(op)
        tableau[p, 2n+1] = bitrand(1)[1]
        outcome = tableau[p,2n+1]
        expVal = 0
    else
        tableau[2n+1, :] .= 0
        @inbounds for i = 1:n
            if check_commute(tableau[i,:], op, n) == false
                tableau = rowsum(tableau, 2n+1, i+n, n)
            end
        end
        outcome = tableau[2n+1, 2n+1]
        expVal = (-1.)^outcome
    end
    return tableau, outcome, expVal
end

function measure_operator_mixed(tableau::BitArray, op::BitArray, n::Int)
    # Return updated tableau, potentially random measurement outcome, and expectation value
    if sum(tableau[1:(2n),2n+2]) == 0
        return measure_operator(tableau, op, n)
    end

    # Case 1
    p = 0
    for i = (n+1):(2n)
        if (tableau[i,2n+2] == 0) && (check_commute(tableau[i,:], op, n) == false)
            p = i
            break
        end
    end
    if p > 0
        for i = 1:2n
            if (i != p) && (check_commute(tableau[i, :], op, n) == false)
                tableau = rowsum(tableau, i, p, n)
            end
        end
        tableau[p - n, :] .= copy(tableau[p, :])
        tableau[p, :] .= op
        tableau[p, 2n+1] = bitrand(1)[1]
        outcome = tableau[p, 2n+1]
        return tableau, outcome, 0
    end

    # Case 2
    p = 0
    for i = 1:(2*n)
        if (tableau[i,2n+2] == 1) && (check_commute(tableau[i, :], op, n) == false)
            p = i
            break
        end
    end
    if p == 0
        tableau[2n+1, :] .= 0
        for i = 1:n
            if (check_commute(tableau[i, :], op, n) == false)
                tableau = rowsum(tableau, 2n+1, i+n, n)
            end
        end
        outcome = tableau[2n+1, 2n+1]
        return tableau, outcome, (-1.)^outcome
    else
        # Case 3
        m = p
        if m > n
            m_bar = m - n
        else
            m_bar = m + n
        end
        for i = 1:(2n)
            if (i != m) && (check_commute(tableau[i,:], op, n) == false)
                tableau = rowsum(tableau, i, m, n)
            end
        end
        tableau[m_bar,:] .= copy(tableau[m,:])
        tableau[m, :] .= copy(op)
        tableau[m, 2n+1] = bitrand(1)[1]
        outcome = tableau[m, 2n+1]
        tableau[m, 2n+2] = 0
        tableau[m_bar, 2n+2] = 0
        if m_bar > m
            swap_row = copy(tableau[m,:])
            tableau[m,:] .= copy(tableau[m_bar, :])
            tableau[m_bar,:] .= swap_row
        end
        return tableau, outcome, 0
    end

end

function measure_zz(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1 + n] = 1
    op[site2 + n] = 1
    tableau, outcome, expVal = measure_operator_mixed(tableau, op, n)
    return tableau, outcome, expVal
end

function measure_x_string(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1:site2] .= 1
    tableau, outcome, expVal = measure_operator_mixed(tableau, op, n)
    return tableau, outcome, expVal
end

function measure_z_string(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[(site1+n):(site2+n)] .= 1
    tableau, outcome, expVal = measure_operator_mixed(tableau, op, n)
    return tableau, outcome, expVal
end

function measure_x(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    tableau, outcome, expVal = measure_operator_mixed(tableau, op, n)
    return tableau, outcome, expVal
end

function measure_z(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site + n] = 1
    tableau, outcome, expVal = measure_operator_mixed(tableau, op, n)
    return tableau, outcome, expVal
end

function measure_y(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    op[site + n] = 1
    tableau, outcome, expVal = measure_operator_mixed(tableau, op, n)
    return tableau, outcome, expVal
end

function initialize_tableau(n::Int)
    # Tableau has 2n+2 columns and 2n+1 rows.
    # First n rows for destabilizers
    # Next n rows for stabilizers
    # Last row for workspace
    # 2n+1 column for phase
    # 2n+2 column for marking when dephasing
    tableau = falses(2n+1, 2n+2)
    for i = 1:(2n)
        tableau[i,i] = 1
    end
    return tableau
end

function hadamard(tableau::BitArray, site::Int, n::Int)
    @inbounds @simd for i = 1:(2*n)
        myX = tableau[i,site]
        myZ = tableau[i,(site+n)]
        myR = tableau[i,2n+1]
        tableau[i,2n+1] = myR ⊻ (myX * myZ) # (myX ⊻ myZ ⊻ 1)
        tableau[i,site] = myZ
        tableau[i,(site+n)] = myX
    end
    return tableau
end

function phase(tableau::BitArray, site::Int, n::Int)
    for i = 1:(2*n)
        tableau[i,2n+1] = tableau[i,2n+1] ⊻ (tableau[i,site] * tableau[i,(site + n)])
        tableau[i,(site + n)] = tableau[i,(site + n)] ⊻ tableau[i,site]
    end
    return tableau
end

function cnot(tableau::BitArray, control_site::Int, target_site::Int, n::Int)
    for i = 1:(2n)
        tableau[i,2n+1] = tableau[i,2n+1] ⊻ (tableau[i, control_site] * tableau[i,control_site+n]) ⊻ (tableau[i,target_site] ⊻ tableau[i, control_site + n] ⊻ 1)
        tableau[i,target_site] = tableau[i,target_site] ⊻ tableau[i,control_site]
        tableau[i,control_site + n] = tableau[i,control_site + n] ⊻ tableau[i,target_site + n]
    end
    return tableau
end

function compute_bipartite_EE(tableau::BitArray, n::Int)
    midpoint = (n ÷ 2)
    my_EE = n - midpoint
    @inbounds @simd for i = 1:n
        if sum(tableau[i+n, (midpoint+1):n] .| tableau[i+n, (midpoint+n+1):2n]) == 0 && tableau[i+n,2n+2]==0
            my_EE -= 1
        end
    end
    return my_EE
end

function dephase_y(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    op[site + n] = 1
    tableau = dephase_operator(tableau, op, n)
    return tableau
end

function dephase_z(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site + n] = 1
    tableau = dephase_operator(tableau, op, n)
    return tableau
end

function dephase_zz(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1 + n] = 1
    op[site2 + n] = 1
    tableau = dephase_operator(tableau, op, n)
    return tableau
end

function dephase_operator(tableau::BitArray, op::BitArray, n::Int)
    p = 0
    for i = (n+1):(2n)
        if tableau[i,2n+2] == 0 && (check_commute(tableau[i,:], op, n) == false)
            p = i
            break
        end
    end

    if p == 0
        return tableau
    end

    for i = 1:(2n)
        if (i != p) && (check_commute(tableau[i,:], op, n) == false)
            tableau = rowsum(tableau, i, p, n)
        end
    end
    tableau[p-n, :] .= tableau[p,:]
    tableau[p, :] .= op

    tableau[p-n, 2n+2] = 1
    tableau[p, 2n+2] = 1

    return tableau

end

function simple_dephase(tableau::BitArray, site::Int, n::Int)
    p = 0
    for i = (n+1):(2n)
        if tableau[i,site+n] == 1 && tableau[i,2n+2] == 0
            p = i
            break
        end
    end
    if p == 0
        return tableau
    end
    for i = 1:(2n)
        if (i != p) && tableau[i,site+n] == 1
            tableau = rowsum(tableau, i, p, n)
        end
    end

    # Copy the anticommuting stabilizer to row p-n and set row p to X_{site}
    tableau[p-n, :] .= tableau[p,:]
    tableau[p, :] .= 0
    tableau[p, site] = 1

    # Mark row as being outside the stabilizer / destabilizer
    tableau[p-n, 2n+2] = 1
    tableau[p, 2n+2] = 1

    # Assertions -- for debugging
    # for i = 1:2n
    #     if (i != p) && (i != p-n)
    #         @assert check_commute(tableau[i,:], tableau[p,:], n) == true
    #         @assert check_commute(tableau[i,:], tableau[p-n,:], n) == true
    #     end
    # end
    # @assert check_commute(tableau[p,:], tableau[p-n,:], n) == false

    return tableau
end

function compute_Z2_rank(M)
    n = size(M)[1]
    nA = size(M)[2]
    i = 1
    for j = 1:nA
        p = findfirst(M[i:n, j])
        if p != nothing
            p = p + i - 1
            tmp = copy(M[i,:])
            M[i, :] .= M[p, :]
            M[p, :] .= tmp
            for k = 1:n# (i+1):n
                if M[k,j] == 1 && i != k
                    M[k,:] .= (M[k,:] .⊻ M[i,:])
                end
            end
            i += 1
        end
    end
    my_rank = 0
    for i = 1:n
        if findfirst(M[i, :]) != nothing
            my_rank += 1
        end
    end
    return my_rank
end

function compute_EE_generic(tableau::BitArray, A::Array{Int64}, n::Int)
    A_extended = union(A, A.+n)
    projected_stabilizers = copy(tableau[(n+1):(2n),A_extended])
    my_rank = compute_Z2_rank(copy(projected_stabilizers))
    my_EE = length(A) - my_rank
    return my_EE
end

function compute_mutual_info(tableau::BitArray, A::Array{Int64}, B::Array{Int64}, n::Int)
    sA = compute_EE_generic(tableau, A, n)
    sB = compute_EE_generic(tableau, B, n)
    sAB = compute_EE_generic(tableau, union(A,B), n)
    return (sA + sB - sAB)
end

function compute_tripartite_info(tableau::BitArray, A::Array{Int64}, B::Array{Int64}, C::Array{Int64}, n::Int)
    sA = compute_EE_generic(tableau, A, n)
    sB = compute_EE_generic(tableau, B, n)
    sC = compute_EE_generic(tableau, C, n)
    sAB = compute_EE_generic(tableau, union(A,B), n)
    sAC = compute_EE_generic(tableau, union(A,C), n)
    sABC = compute_EE_generic(tableau, union(A,B,C), n)
    return (sA + sB + sC - sAB - sAC + sABC)
end

println("loaded 2")
