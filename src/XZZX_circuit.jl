"""
    site_to_op(site, n)

Construct the XZZX operator starting at the specified site for a system of n
    qubits.
"""
function site_to_op(site::Int, n::Int)
    op = falses(2n+2)
    cluster = [mod(i-1, n) + 1 for i=site:(site+5)]
    op[cluster[1]] = 1
    op[cluster[2] + n] = 1
    op[cluster[3] + n] = 1
    op[cluster[4]] = 1
    return op
end

"""
    site_to_SG_op(site1, site2, n)

Construct the spin glass operator M_{site1} M_{site1+1} ... M_{site2} where
M = XZZX and there are n qubits in total. We assume that site1 and site2 are
sufficiently well separated that there is no overlap between the XYX at either
end, which will be true for the typical use case.
"""
function sites_to_SG_op(site1::Int, site2::Int, n::Int)
    op = falses(2n+2)
    op[site1] = 1
    op[site2] = 1
    op[mod(site1 - 2, n) + 1] = 1
    op[mod(site1, n) + 1] = 1
    op[mod(site2 - 2, n) + 1] = 1
    op[mod(site2, n) + 1] = 1
    op[site1 + n] = 1
    op[site2 + n] = 1
    return op
end

"""
    measure_spin_glass(t, n)

For an n-qubit stabilizer state given by the tableau, compute the spin glass
order parameter for the stabilizers (M = XZZX). See e.g. arXiv:2108.04274 for
discussion of the order parameter.
"""
function measure_spin_glass(tableau::BitArray, n::Int)
    A = [i for i=1:(n÷8)] .+ (n÷4)
    B = A .+ (n÷2)
    my_matrix = zeros(length(A), length(B))
    for i in A
        for j in B
            expect_val = operator_expectation_mixed(copy(tableau), sites_to_SG_op(i,j,n), n)[3]
            my_matrix[i - A[1] + 1,j - B[1] + 1] = expect_val
        end
    end
    return my_matrix
end

"""
    update_circuit_stabilizers(t, stabs, n, pStab, PBC=PBC)

For an n-qubit stabilizer state given by the tableau, update the state for a
layer of the circuit involving only the stabilizers M_i (XZZX). The stabilizers
are given in an array stabs. Each stabilizer is measured with probability pStab.
Optional argument PBC determines the boundaryt conditions.
"""
function update_circuit_stabilizers!(tableau::BitArray, stabs, n::Int, pStab::Float64; PBC=true)
    for j = 1:n
        if ((j + 3) > n) && !PBC
            break
        end
        if rand() < pStab
            measure_operator_mixed!(tableau, stabs[j], n)
        end
    end
end

"""
    update_circuit_errors(t, n, pX, pY, pZ)

For an n-qubit stabilizer state given by the tableau, update the state for a
layer of the circuit involving only the errors (X, Y, Z). The measurement
probabilities are given by pX, pY, and pZ.
"""
function update_circuit_errors!(tableau::BitArray, n::Int, pX::Float64, pY::Float64, pZ::Float64)
    for j = 1:n
        r = rand()
        if r < pX
            measure_x!(tableau, j, n)
        elseif r < (pX + pY)
            measure_y!(tableau, j, n)
        elseif r < (pX + pY + pZ)
            measure_z!(tableau, j, n)
        end
    end
end

"""
    update_circuit_single_layer(t, n, pX, pY, pZ, stabs, PBC=PBC)

For an n-qubit stabilizer state given by the tableau, update the state for a
layer of the circuit involving only the errors (X, Y, Z). The measurement
probabilities are given by pX, pY, and pZ. Stabilizers M=XZZX are stored in the
list "stabs" and are measured with probability 1 - pX - pY - pZ. Optional
argument PBC determines the boundary conditions applied for stabilizer
measurements.

Unlike the other update_circuit functions, here we adopt the update scheme from
arXiv:2108.04274, wherein a single timestep consists of n randomly positioned
measurements, with the appropriate probability for the measured operator.
"""
function update_circuit_single_layer!(tableau::BitArray, n::Int, pX::Float64, pY::Float64, pZ::Float64, stabs; PBC=true)
    i = 0
    while i < n
        site = rand([k for k = 1:n])
        r = rand()
        if r < (pX)
            measure_x!(tableau, site, n)
        elseif r < (pX + pY)
            measure_y!(tableau, site, n)
        elseif r < (pX + pY + pZ)
            measure_z!(tableau, site, n)
        else
            if (site + 3 < n) || PBC
                measure_operator_mixed!(tableau, stabs[site], n)
            else
                i -= 1
            end
        end
        i += 1
    end
end

"""
    run_circuit(n, pX, pY, pZ, sat_steps, avg_steps, avg_interval, PBC=PBC)

Run a simulation of an n-qubit measurement-only XZZX circuit. Error measurement
probabilities are given by pX, pY, pZ. The circuit is run for sat_step layers
to reach the steady state. In the steady state, avg_step layers are averaged
over, with avg_interval layers between each value recorded. Boundary conditions
are set by optional argument PBC.

Circuit structure is alternating between stabilizer layers and error layers.

Returns the average of the following entanglement measures:
    1) Half-chain entanglement entropy
    2) Antipodal mutual information
    3) S_topo^c
    4) S_topo^t
    5) S_topo^q
"""
function run_circuit(n::Int, pX::Float64, pY::Float64, pZ::Float64, sat_steps::Int, avg_steps::Int, avg_interval::Int; PBC=true)
    pStab = 1 - pX - pY - pZ

    stab_list = [site_to_op(i, n) for i = 1:n]

    my_tableau = initialize_tableau(n)

    # Alternating layer circuit evolution to reach the steady state
    step_num = 0
    for i=1:sat_steps
        if step_num % 2 == 0
            update_circuit_stabilizers!(my_tableau, stab_list, n, pStab, PBC=PBC)
        else
            update_circuit_errors!(my_tableau, n, pX, pY, pZ)
        end
        step_num += 1
    end

    # Compute the average entanglement measures in the steady state
    data_history = zeros(avg_steps, 5)
    for i = 1:avg_steps
        data_history[i,1] = compute_bipartite_EE(copy(my_tableau), n)
        data_history[i,2] = compute_antipodal_mutual(copy(my_tableau), n)
        data_history[i,3] = compute_Stopo_c(copy(my_tableau), n)
        data_history[i,4] = compute_Stopo_t(copy(my_tableau), n)
        data_history[i,5] = compute_Stopo_q(copy(my_tableau), n)
        for j = 1:avg_interval
            if step_num % 2 == 0
                update_circuit_stabilizers!(my_tableau, stab_list, n, pStab, PBC=PBC)
            else
                update_circuit_errors!(my_tableau, n, pX, pY, pZ)
            end
            step_num += 1
        end
    end

    mean_data = mean(data_history,dims=1)[1,:]

    return mean_data
end

"""
    run_circuit_single_layer(n, pX, pY, pZ, sat_steps, avg_steps, avg_interval, PBC=PBC)

Run a simulation of an n-qubit measurement-only XZZX circuit. Error measurement
probabilities are given by pX, pY, pZ. The circuit is run for sat_step layers
to reach the steady state. In the steady state, avg_step layers are averaged
over, with avg_interval layers between each value recorded. Boundary conditions
are set by optional argument PBC.

Unlike run_circuit, this update scheme does all measurements in a single layer,
as in arXiv:2004.07243

    Returns the average of the following entanglement measures:
        1) Half-chain entanglement entropy
        2) Antipodal mutual information
        3) S_topo^c
        4) S_topo^t
        5) S_topo^q
"""
function run_circuit_single_layer(n::Int, pX::Float64, pY::Float64, pZ::Float64, sat_steps::Int, avg_steps::Int, avg_interval::Int; PBC=true)
    pStab = 1. - pX - pY - pZ

    stab_list = [site_to_op(i, n) for i = 1:n]

    my_tableau = initialize_tableau(n)

    # Reach the steady state
    for i=1:(sat_steps)
        update_circuit_single_layer!(my_tableau, n, pX, pY, pZ, stab_list, PBC=PBC)
    end

    # Measure the entanglement metrics in the steady state
    data_history = zeros(avg_steps, 5)
    for i = 1:avg_steps
        data_history[i,1] = compute_bipartite_EE(copy(my_tableau), n)
        data_history[i,2] = compute_antipodal_mutual(copy(my_tableau), n)
        data_history[i,3] = compute_Stopo_c(copy(my_tableau), n)
        data_history[i,4] = compute_Stopo_t(copy(my_tableau), n)
        data_history[i,5] = compute_Stopo_q(copy(my_tableau), n)
        for j = 1:(avg_interval)
            update_circuit_single_layer!(my_tableau, n, pX, pY, pZ, stab_list, PBC=PBC)
        end
    end

    mean_data = mean(data_history,dims=1)[1,:]

    return mean_data
end

"""
    run_circuit_EE_growth(n, pX, pY, pZ, sat_steps, avg_steps, avg_interval, PBC=PBC)

Run a simulation of an n-qubit measurement-only XZZX circuit. Error measurement
probabilities are given by pX, pY, pZ. Boundary conditions are set by optional
argument PBC.

The time-dependence of the half-chain entanglement entropy is recorded for
sat_step circuit layers.
"""
function run_circuit_EE_growth(n::Int, pX::Float64, pY::Float64, pZ::Float64, sat_steps::Int; PBC=true)
    pStab = 1. - pX - pY - pZ

    stab_list = [site_to_op(i, n) for i = 1:n]

    my_tableau = initialize_tableau(n)

    ee_history = zeros(sat_steps)

    for i = 1:sat_steps
        ee_history[i] = compute_bipartite_EE(copy(my_tableau), n)
        if i % 2 == 0
            update_circuit_stabilizers!(my_tableau, stab_list, n, pStab, PBC=PBC)
        else
            update_circuit_errors!(my_tableau, n, pX, pY, pZ)
        end
    end

    return ee_history
end

"""
    run_circuit_EE_scaling(n, pX, pY, pZ, sat_steps, avg_steps, avg_interval, PBC=PBC)

Run a simulation of an n-qubit measurement-only XZZX circuit. Error measurement
probabilities are given by pX, pY, pZ. The circuit is run for sat_step layers
to reach the steady state. In the steady state, avg_step layers are averaged
over, with avg_interval layers between each value recorded. Boundary conditions
are set by optional argument PBC.

Circuit structure is alternating between stabilizer layers and error layers.

Returns the spatial scaling of the entanglement entropy averaged in the steady
state.
"""
function run_circuit_EE_scaling(n::Int, pX::Float64, pY::Float64, pZ::Float64, sat_steps::Int, avg_steps::Int, avg_interval::Int; PBC=true)
    pStab = 1 - pX - pY - pZ

    stab_list = [site_to_op(i, n) for i = 1:n]

    my_tableau = initialize_tableau(n)

    # Reach the steady state
    step_num = 0
    for j = 1:(sat_steps*n)
        if step_num % 2 == 0
            update_circuit_stabilizers!(my_tableau, stab_list, n, pStab, PBC=PBC)
        else
            update_circuit_errors!(my_tableau, n, pX, pY, pZ)
        end
        step_num += 1
    end

    # Average the entanglement scaling over the steady state.
    ee_scaling_history = zeros(avg_steps, n)
    for i = 1:avg_steps
        ee_scaling_history[i,:] = compute_entanglement_scaling(copy(my_tableau), n)
        for j = 1:(avg_interval)
            if step_num % 2 == 0
                update_circuit_stabilizers!(my_tableau, stab_list, n, pStab, PBC=PBC)
            else
                update_circuit_errors!(my_tableau, n, pX, pY, pZ)
            end
            step_num += 1
        end
    end

    mean_ee_scaling = mean(ee_scaling_history, dims=1)[1,:]

    return mean_ee_scaling

end

"""
    purify_mixed_state(n, pX, pY, pZ, Tmax, PBC=PBC)

Run a simulation for Tmax layers on an n-qubit measurement-only XZZX circuit.
Starting from an initial maximally mixed state, the purity and residual entropy
are recorded at every step. Probabilities pX, pY, pZ and boundary conditions PBC
control the update for the circuit.

Circuit structure is alternating between stabilizer layers and error layers.
"""
function purify_mixed_state(n::Int, pX::Float64, pY::Float64, pZ::Float64, Tmax::Int; PBC=false)
    pStab = 1. - pX - pY - pZ

    stab_list = [site_to_op(i, n) for i = 1:n]

    my_tableau = initialize_tableau(n)

    # Prepare a maximally mixed initial state
    for i = 1:n
        dephase_x!(my_tableau, i, n)
    end

    # Evolve the circuiit and save the residual entropy density and purity.
    data_history = zeros(Tmax, 2)
    data_history[:,1] .= 1. # Default for when break condition is reached
    for i = 1:Tmax
        data_history[i,2] = compute_residual_entropy(my_tableau, n)
        data_history[i,1] = 2 .^ (-data_history[i,2])
        if i % 2 == 1
            update_circuit_stabilizers!(my_tableau, stab_list, n, pStab, PBC=PBC)
        else
            update_circuit_errors!(my_tableau, n, pX, pY, pZ)
        end
        # Break if the state is fully purified
        if data_history[i,2] == 0
            break
        end
    end

    return data_history
end
