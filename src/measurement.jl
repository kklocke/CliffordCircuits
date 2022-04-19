
"""
    measure_operator!(t, op, n)

Do an (in-place) measurement of an n-qubit operator op on a stabilizer specified
by the tableau. The tableau should be for a pure state (all entries in the last
column must be zero).

Crucially, the expectation value is *not* returned.
"""
function measure_operator!(tableau::BitArray, op::BitArray, n::Int)
    p = 0
    @inbounds for i = (n+1):(2n)
        if !check_commute(tableau[i,:], op, n)
            p = i
            break
        end
    end

    if p > 0
        @inbounds for i = 1:2n
            if (i != p) && !check_commute(tableau[i,:], op, n)
                tableau = rowsum(tableau, i, p, n)
            end
        end
        tableau[p-n, :] .= tableau[p,:]
        tableau[p, :] .= op
        tableau[p, 2n+1] = bitrand(1)[1]
    end
end

"""
    measure_operator_mixed!(t, op, n)

Do an in-place measurement of an n-qubit operator op on a stabilizer specified
by the tableau. The tableau may be for a mixed state.

The expectation value is *not* returned.
"""
function measure_operator_mixed!(tableau::BitArray, op::BitArray, n::Int)
    # Return updated tableau, potentially random measurement outcome, and expectation value
    if sum(tableau[1:(2n),2n+2]) == 0
        measure_operator!(tableau, op, n)
        return
    end

    # Case 1
    p = 0
    for i = (n+1):(2n)
        if (tableau[i,2n+2] == 0) && !check_commute(tableau[i,:], op, n)
            p = i
            break
        end
    end
    if p > 0
        for i = 1:2n
            if (i != p) && !check_commute(tableau[i, :], op, n)
                tableau = rowsum(tableau, i, p, n)
            end
        end
        tableau[p - n, :] .= tableau[p, :]
        tableau[p, :] .= op
        tableau[p, 2n+1] = bitrand(1)[1]
        return
    end

    # Case 2
    p = 0
    for i = 1:(2*n)
        if (tableau[i,2n+2] == 1) && !check_commute(tableau[i, :], op, n)
            p = i
            break
        end
    end
    if p > 0
        # Case 3
        m = p
        if m > n
            m_bar = m - n
        else
            m_bar = m + n
        end
        for i = 1:(2n)
            if (i != m) && !check_commute(tableau[i,:], op, n)
                tableau = rowsum(tableau, i, m, n)
            end
        end
        tableau[m_bar,:] .= tableau[m,:]
        tableau[m, :] .= op
        tableau[m, 2n+1] = bitrand(1)[1]
        tableau[m, 2n+2] = 0
        tableau[m_bar, 2n+2] = 0
        if m_bar > m
            rowswap!(tableau, m, m_bar)
        end
    end
end

"""
    operator_expectation(t, op, n)

Compute the expectation value of an n-qubit operator op given a pure stablizer
state specified by the tableau. This function does *not* act in place, so that
you may pass a copy of the tableau to get the expectation value without altering
the state.
"""
function operator_expectation(tableau::BitArray, op::BitArray, n::Int)
    p = 0
    @inbounds for i = (n+1):(2n)
        if !check_commute(tableau[i,:], op, n)
            p = i
            break
        end
    end

    if p > 0
        @inbounds for i = 1:2n
            if (i != p) && !check_commute(tableau[i,:], op, n)
                tableau = rowsum(tableau, i, p, n)
            end
        end
        tableau[p-n, :] .= tableau[p,:]
        tableau[p, :] .= op
        tableau[p, 2n+1] = bitrand(1)[1]
        outcome = tableau[p,2n+1]
        expVal = 0
    else
        tableau[2n+1, :] .= 0
        @inbounds for i = 1:n
            if !check_commute(tableau[i,:], op, n)
                tableau = rowsum(tableau, 2n+1, i+n, n)
            end
        end
        outcome = tableau[2n+1, 2n+1]
        expVal = (-1.)^outcome
    end
    return tableau, outcome, expVal
end

"""
    operator_expectation_mixed(t, op, n)

Compute the expectation value of an n-qubit operator op given a mixed stablizer
state specified by the tableau. This function does *not* act in place, so that
you may pass a copy of the tableau to get the expectation value without altering
the state.

"""
function operator_expectation_mixed(tableau::BitArray, op::BitArray, n::Int)
    # Return updated tableau, potentially random measurement outcome, and expectation value
    if sum(tableau[1:(2n),2n+2]) == 0
        return operator_expectation(tableau, op, n)
    end

    # Case 1
    p = 0
    for i = (n+1):(2n)
        if (tableau[i,2n+2] == 0) && !check_commute(tableau[i,:], op, n)
            p = i
            break
        end
    end
    if p > 0
        for i = 1:2n
            if (i != p) && !check_commute(tableau[i, :], op, n)
                tableau = rowsum(tableau, i, p, n)
            end
        end
        tableau[p - n, :] .= tableau[p, :]
        tableau[p, :] .= op
        tableau[p, 2n+1] = bitrand(1)[1]
        outcome = tableau[p, 2n+1]
        return tableau, outcome, 0
    end

    # Case 2
    p = 0
    for i = 1:(2*n)
        if (tableau[i,2n+2] == 1) && !check_commute(tableau[i, :], op, n)
            p = i
            break
        end
    end
    if p == 0
        tableau[2n+1, :] .= 0
        for i = 1:n
            if !check_commute(tableau[i, :], op, n)
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
            if (i != m) && !check_commute(tableau[i,:], op, n)
                tableau = rowsum(tableau, i, m, n)
            end
        end
        tableau[m_bar,:] .= tableau[m,:]
        tableau[m, :] .= op
        tableau[m, 2n+1] = bitrand(1)[1]
        outcome = tableau[m, 2n+1]
        tableau[m, 2n+2] = 0
        tableau[m_bar, 2n+2] = 0
        if m_bar > m
            !rowswap(tableau, m, m_bar)
        end
        return tableau, outcome, 0
    end
end

"""
    measure_zz(t, site1, site2, n)

Measure in-place the operator Z_{site 1} Z_{site 2} on an n-qubit stabilizer
state given by the tableau.
"""
function measure_zz!(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1 + n] = 1
    op[site2 + n] = 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    measure_xx(t, site1, site2, n)

Measure in-place the operator X_{site 1} X_{site 2} on an n-qubit stabilizer
state given by the tableau.
"""
function measure_xx!(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1] = 1
    op[site2] = 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    measure_x_string(t, site1, site2, n)

This function generalizes measure_xx by allowing for in-place measurement of a
string  X_{site1} X_{site1+1} ... X_{site2-1} X_{site2} on an n-qubit stabilizer
state given by the tableau.
"""
function measure_x_string!(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1:site2] .= 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    measure_z_string(t, site1, site2, n)

This function generalizes measure_zz by allowing for in-place measurement of a
string  Z_{site1} Z_{site1+1} ... Z_{site2-1} Z_{site2} on an n-qubit stabilizer
state given by the tableau.
"""
function measure_z_string!(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[(site1+n):(site2+n)] .= 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    measure_x(t, site, n)

Measure in-place the operator X_{site} for an n-qubit stabilizer state given by
the tableau.
"""
function measure_x!(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    measure_z(t, site, n)

Measure in-place the operator Z_{site} for an n-qubit stabilizer state given by
the tableau.
"""
function measure_z!(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site + n] = 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    measure_y(t, site, n)

Measure in-place the operator Y_{site} for an n-qubit stabilizer state given by
the tableau.
"""
function measure_y!(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    op[site + n] = 1
    op[2n+1] = 1
    measure_operator_mixed!(tableau, op, n)
end

"""
    dephase_operator!(t, op, n)

For an n-qubit stabilizer state given by the tableau, dephase the operator op.
"""
function dephase_operator!(tableau::BitArray, op::BitArray, n::Int)
    p = 0
    for i = (n+1):(2n)
        if tableau[i,2n+2] == 0 && !check_commute(tableau[i,:], op, n)
            p = i
            break
        end
    end

    if p == 0
        return
    end

    for i = 1:(2n)
        if (i != p) && !check_commute(tableau[i,:], op, n)
            tableau = rowsum(tableau, i, p, n)
        end
    end
    tableau[p-n, :] .= tableau[p,:]
    tableau[p, :] .= op

    tableau[p-n, 2n+2] = 1
    tableau[p, 2n+2] = 1
end

"""
    dephase_x!(t, op, n)

For an n-qubit stabilizer state given by the tableau, dephase operator X_{site}.
"""
function dephase_x!(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    dephase_operator!(tableau, op, n)
end

"""
    dephase_y!(t, op, n)

For an n-qubit stabilizer state given by the tableau, dephase operator Y_{site}.
"""
function dephase_y!(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site] = 1
    op[site + n] = 1
    dephase_operator!(tableau, op, n)
end

"""
    dephase_z!(t, op, n)

For an n-qubit stabilizer state given by the tableau, dephase operator Z_{site}.
"""
function dephase_z!(tableau::BitArray, site::Int, n::Int)
    op = falses(2n+2)
    op[site + n] = 1
    dephase_operator!(tableau, op, n)
end

"""
    dephase_zz(t, op, n)

For an n-qubit stabilizer state given by the tableau, dephase operator
Z_{site1}, Z_{site2}.
"""
function dephase_zz!(tableau::BitArray, site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1 + n] = 1
    op[site2 + n] = 1
    dephase_operator!(tableau, op, n)
end
