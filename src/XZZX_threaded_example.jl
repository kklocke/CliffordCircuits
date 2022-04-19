@everywhere include("CliffordCircuits.jl")
@everywhere include("XZZX_circuit.jl")
using Distributed

# Set the number of qubits n
global n = 18
if length(ARGS) >= 1
    n = parse(Int64, ARGS[1])
end
@show n

# Set qX
global qX = 0.25
if length(ARGS) >= 2
    qX = parse(Float64, ARGS[2])
end
@show qX

# Set qY
global qY = 0.25
if length(ARGS) >= 3
    qY = parse(Float64, ARGS[3])
end
@show qY

qZ = 1. - qX - qY

# Number of repetitions to average over for each p_M^{Stab}
global num_reps = 2
if length(ARGS) >= 4
    num_reps = parse(Int64, ARGS[4])
end
@show num_reps

# Default values for reaching the steady state and averaging.
sat_steps = max(200,5*n)
avg_interval = 31
avg_steps = 100

# Desired range for p_M^{Stab}
p_range = [0., 0.1, 0.2, 0.3, 0.35, 0.4, 0.425, 0.45, 0.475, 0.4875, 0.49375, 0.5, 0.50625, 0.5125, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0]
p_range_pmap = [p for p in p_range for j=1:num_reps]

myRes = pmap(p->run_circuit(n, qX*(1-p), qY*(1-p), qZ*(1-p), sat_steps, avg_steps, avg_interval, PBC=false), p_range_pmap)

all_dat = [elem for elem in myRes]
all_dat = reshape(cat(all_dat...,dims=2),5,num_reps,length(p_range))
all_dat = permutedims(all_dat,[3,2,1])

# Compute the mean and standard deviation for each value of p_M^{Stab}
m_dat = mean(all_dat,dims=2)[:,1,:]
s_dat = std(all_dat,dims=2)[:,1,:]

# Save p_M^{Stab}, the mean, and the standard deviation
fname = string("Results/code513_qX", round(qX,digits=3), "_qY",round(qY,digits=3),"_N", n, "_data_OBC.txt")
to_write = cat(p_range, m_dat, s_dat, dims=2)
open(fname, "w") do io
    writedlm(io, to_write, ',')
end
