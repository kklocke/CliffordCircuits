"""
    compute_Z2_rank(M)

For a matrix M (not necessarily square), compute the rank of the matrix with
respect to the field Z_2. Calculation proceeds via Gaussian elimination.
"""
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
            for k = 1:n
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

"""
    compute_EE_generic(t, A, n)

For an n-qubit stabilizer state given by the tableau, compute the (von Neumann)
entanglement entropy of the region A.

Argument A should be a list of integers specifying the sites contained in A.
"""
function compute_EE_generic(tableau::BitArray, A::Array{Int64}, n::Int)
    A_extended = union(A, A.+n)
    mask = [i for i = (n+1):(2n)  if tableau[i, 2n+2] == 0]
    projected_stabilizers = copy(tableau[mask, A_extended])
    my_rank = compute_Z2_rank(copy(projected_stabilizers))
    my_EE = my_rank - length(A)
    return my_EE
end

"""
    compute_bipartite_EE(t, n)

For an n-qubit stabilizer state given by the tableau, compute the entanglement
entropy between the left and right half of the system.
"""
function compute_bipartite_EE(tableau::BitArray, n::Int)
    midpoint = (n ÷ 2)
    A = [i for i = 1:midpoint]
    return compute_EE_generic(tableau, A, n)
end

"""
    compute_entanglement_scaling(t, n)

For an n-qubit stabilizer state given by the tableau, compute the entanglement
scaling by calculating the entanglement entropy of all intervals [1, ... , l].
"""
function compute_entanglement_scaling(tableau::BitArray, n::Int)
    res = zeros(n)
    stabs = copy(tableau[(n+1):(2n), :])
    for i=1:(n-1)
        A = [j for j = 1:i]
        res[i] = compute_EE_generic(copy(tableau), A, n)
    end
    return res
end

"""
    compute_mutual_info(t, A, B, n)

For an n-qubit stabilizer state given by the tableau, compute the mutual
information ``I_2(A,B) = S_A + S_B - S_{AB}`` between subsystem A and B.
"""
function compute_mutual_info(tableau::BitArray, A::Array{Int64}, B::Array{Int64}, n::Int)
    sA = compute_EE_generic(tableau, A, n)
    sB = compute_EE_generic(tableau, B, n)
    sAB = compute_EE_generic(tableau, union(A,B), n)
    return (sA + sB - sAB)
end

"""
    compute_antipodal_mutual(t, n)

For an n-qubit stabilizer state given by the tableau, compute the mutual
information between regions A and B where A consists of the first L/8 sites,
and B is displaced from A by L/2.
"""
function compute_antipodal_mutual(tableau::BitArray, n::Int)
    A = [i for i = 1:(n ÷ 8)]
    B = A .+ (n ÷ 2)
    return compute_mutual_info(tableau, A, B, n)
end

"""
    compute_Stopo(t, A, B, C, n)

For an n-qubit stabilizer state given by the tableau, compute the topological
entanglement entropy with respect to subsystems A, B, C.
In particular we take ``S_{topo} = S_{AB} + S_{BC} - S_B - S{ABC}``.
"""
function compute_Stopo(tableau::BitArray, A::Array{Int64}, B::Array{Int64}, C::Array{Int64}, n::Int)
    sA = compute_EE_generic(tableau, A, n)
    sB = compute_EE_generic(tableau, B, n)
    sC = compute_EE_generic(tableau, C, n)
    sAB = compute_EE_generic(tableau, union(A,B), n)
    sAC = compute_EE_generic(tableau, union(A,C), n)
    sBC = compute_EE_generic(tableau, union(B,C), n)
    sABC = compute_EE_generic(tableau, union(A,B,C), n)
    return (sAB + sBC - sB - sABC)
end


"""
    compute_Stopo_c(t,n)

For an n-qubit stabilizer state given by the tableau, compute S_topo^c to
identify the critical point and log-law coefficient. This computes S_topo on
a partitioning of the system into four equally sized regions A, B, C, D.
"""
function compute_Stopo_c(tableau::BitArray, n::Int)
    A = [i for i = 1:(n ÷ 4)]
    B = A .+ (n ÷ 4)
    C = B .+ (n ÷ 4)
    return compute_Stopo(tableau, A, B, C, n)
end

"""
    compute_Stopo_t(t,n)

For an n-qubit stabilizer state given by the tableau, compute S_topo^t. This
computes S_topo on a partitioning of the system into three equally sized regions
A, B, C.
"""
function compute_Stopo_t(tableau::BitArray, n::Int)
    A = [i for i = 1:(n ÷ 3)]
    B = A .+ (n ÷ 3)
    C = setdiff([i for i=1:n], union(A,B))
    return compute_Stopo(tableau, A, B, C, n)
end

"""
    compute_Stopo_q(t,n)

For an n-qubit stabilizer state given by the tableau, compute S_topo^q. This
computes S_topo on a partitioning of the system into four equally sized regions
A, B, D, C.
"""
function compute_Stopo_q(tableau::BitArray, n::Int)
    A = [i for i = 1:(n ÷ 4)]
    B = A .+ (n ÷ 4)
    C = B .+ (n ÷ 2)
    return compute_Stopo(tableau, A, B, C, n)
end

"""
    compute_residual_entropy(t, n)

For an n-qubit stabilizer state given by the tableau, compute the residual
entropy of the state. This is just the number of generators of the mixed space,
which can be found by summing the 2n+2 column of the tableau.
For a pure state the result will be zero.
"""
function compute_residual_entropy(tableau::BitArray, n::Int)
    purity = sum(tableau[(n+1):2n, 2n+2])
    return purity
end
