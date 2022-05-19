"""
    hadamard!(t, site, n)

For an n-qubit stabilizer state given by the tableau, apply the Hadamard gate
to the specified site.
"""
function hadamard!(tableau::BitArray, site::Int, n::Int)
    @inbounds @simd for i = 1:(2*n)
        myX = tableau[i,site]
        myZ = tableau[i,(site+n)]
        myR = tableau[i,2n+1]
        tableau[i,2n+1] = myR ⊻ (myX * myZ) # (myX ⊻ myZ ⊻ 1)
        tableau[i,site] = myZ
        tableau[i,(site+n)] = myX
    end
end

"""
    phase!(t, site, n)

For an n-qubit stabilizer state given by the tableau, apply the phase gate
to the specified site.
"""
function phase!(tableau::BitArray, site::Int, n::Int)
    for i = 1:(2*n)
        tableau[i,2n+1] = tableau[i,2n+1] ⊻ (tableau[i,site] * tableau[i,(site + n)])
        tableau[i,(site + n)] = tableau[i,(site + n)] ⊻ tableau[i,site]
    end
end

"""
    cnot!(t, control_site, target_site, n)

For an n-qubit stabilizer state given by the tableau, apply the CNOT gate, with
the control site and target site as specified by the input.
"""
function cnot!(tableau::BitArray, control_site::Int, target_site::Int, n::Int)
    for i = 1:(2n)
        tableau[i,2n+1] = tableau[i,2n+1] ⊻ (tableau[i, control_site] * tableau[i,control_site+n]) ⊻ (tableau[i,target_site] ⊻ tableau[i, control_site + n] ⊻ 1)
        tableau[i,target_site] = tableau[i,target_site] ⊻ tableau[i,control_site]
        tableau[i,control_site + n] = tableau[i,control_site + n] ⊻ tableau[i,target_site + n]
    end
end

"""
    Z_gate!(t, site, n)

For an n-qubit stabilizer state given by the tableau, apply the Z gate at the
specified site.
"""
function Z_gate!(tableau::BitArray, site::Int, n::Int)
    phase!(tableau, site, n)
    phase!(tableau, site, n)
end

"""
    X_gate!(t, site, n)

For an n-qubit stabilizer state given by the tableau, apply the X gate at the
specified site.
"""
function X_gate!(tableau::BitArray, site::Int, n::Int)
    hadamard!(tableau, site, n)
    Z_gate!(tableau, site, n)
    hadamard!(tableau, site, n)
end

"""
    CZ_gate!(t, control_site, target_site, n)

For an n-qubit stabilizer state given by the tableau, apply the CZ gate, with
the control site and target site as specified by the input.
"""
function CZ_gate!(tableau::BitArray, control_site::Int, target_site::Int, n::Int)
    hadamard!(tableau, target_site, n)
    cnot!(tableau, control_site, target_site, n)
    hadamard!(tableau, target_site, n)
end
