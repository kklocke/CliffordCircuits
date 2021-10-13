include("stabilizer_helper.jl")

function site_to_op_test(n)
    op = falses(2n+2)
    # cluster = [mod(i-1, n) + 1 for i=site:(site+5)]
    # op[cluster[1]] = 1
    # op[cluster[2] + n] = 1
    # op[cluster[3] + n] = 1
    # op[cluster[4]] = 1
    return op
end

function test_func(x)
    return x^2 + test_func_2(x)
end



function site_to_op(site::Int, n::Int)
    op = falses(2n+2)
    cluster = [mod(i-1, n) + 1 for i=site:(site+5)]
    op[cluster[1]] = 1
    op[cluster[2] + n] = 1
    op[cluster[3] + n] = 1
    op[cluster[4]] = 1
    return op
end

function sites_to_SG_op(site1::Int, site2::Int, n::Int)
    # Assuming by default that there is no overlap
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

function measure_spin_glass(tableau::BitArray, n::Int)
    A = [i for i=1:(n÷8)] .+ (n÷4)
    B = A .+ (n÷2)
    my_matrix = zeros(length(A), length(B))
    for i in A
        for j in B
            t, mVal, ez = measure_operator_mixed(copy(tableau), sites_to_SG_op(i,j,n), n)
            # t,z12,ez = measure_zz(copy(tableau), i, j, n)
            my_matrix[i - A[1] + 1,j - B[1] + 1] = ez
        end
    end
    return my_matrix
end

function measure_paramagnetic_X(tableau::BitArray, n::Int)
    A = [i for i=1:(n÷8)] .+ (n÷4)
    B = A .+ (n÷2)
    my_matrix = zeros(length(A), length(B))
    for i in A
        for j in B
            t,x,ex = measure_x_string(copy(tableau), i, j, n)
            my_matrix[i - A[1] + 1, j - B[1] + 1] = ex
        end
    end
    return my_matrix
end

function measure_paramagnetic_Z(tableau::BitArray, n::Int)
    A = [i for i=1:(n÷8)] .+ (n÷4)
    B = A .+ (n÷2)
    my_matrix = zeros(length(A), length(B))
    for i in A
        for j in B
            t, x, ex = measure_z_string(copy(tableau), i, j, n)
            my_matrix[i - A[1] + 1, j - B[1] + 1] = ex
        end
    end
    return my_matrix
end

function measure_mean_Z(tableau::BitArray, n::Int)
    my_mean = 0.
    for i = 1:n
        my_mean += abs(measure_z(copy(tableau), i, n)[3])
    end
    return (my_mean / n)
end

function measure_mean_X(tableau::BitArray, n::Int)
    my_mean = 0.
    for i = 1:n
        my_mean += abs(measure_x(copy(tableau), i, n)[3])
    end
    return (my_mean / n)
end

function measure_entanglement_scaling(tableau::BitArray, n::Int)
    res = zeros(n)
    stabs = copy(tableau[(n+1):(2n), :])
    for i = 1:(n-1)
        for offset = 0:(n-1)
            A = [mod(j + offset - 1, n)+1 for j=1:i]
            myEE = compute_EE_generic(copy(tableau), A, n)
            B = setdiff([i for i=1:n], A)
            myEE2 = compute_EE_generic(copy(tableau), B, n)
            if (myEE2 != myEE)
                @assert myEE == myEE2
            end
            res[i] += (myEE / n)
        end
    end
    return res
end

function measure_antipodal_mutual(tableau::BitArray, n::Int)
    # A = [i for i = 1:(n ÷ 8)]
    A = [1]
    B = A .+ (n ÷ 2)
    return compute_mutual_info(tableau, A, B, n)
end

function measure_tripartite(tableau::BitArray, n::Int)
    A = [i for i = 1:(n ÷ 4)]
    B = A .+ (n ÷ 4)
    C = B .+ (n ÷ 4)
    return compute_tripartite_info(tableau, A, B, C, n)
end

function measure_stabilizer_scaling(tableau::BitArray, n::Int)
    res = zeros(n ÷ 2)
    for i = 4:length(res)
        res[i] = abs(measure_operator(tableau, sites_to_SG_op(1,i,n), n)[3])
    end
    return res
end

function measure_mutual_scaling(tableau::BitArray, n::Int)
    res = zeros(n)
    counts = zeros(n)
    for i = 1:n
        for d = 1:(n-1)
            j = mod(i - 1 + d, n) + 1
            counts[d] += 1
            res[d] += compute_mutual_info(tableau, [i], [j], n)
        end
    end
    for i=1:n
        if res[i] > 0
            res[i] /= counts[i]
        end
    end
    return res
end

function update_circuit(tableau::BitArray, n::Int, pStab::Float64, pX::Float64, pE::Float64, i::Int)
    j_range = 1:n
    j_range = shuffle(j_range)
    if i % 2 == 0
        for j in j_range
            if rand() < pStab
                tableau, mVal, expVal = measure_operator_mixed(tableau, site_to_op(j,n), n)
                # tableau, mVal, expVal = measure_operator_mixed(tableau, stab_list[j], n)
            end
        end
    else
        for j in j_range
            r = rand()
            if r < pX
                tableau, xVal, expVal = measure_x(tableau, j, n)
            elseif r < (pX + pE)
                tableau, zVal, expVal = measure_z(tableau, j, n)
            end
        end
    end
    return tableau
end

function run_circuit(n::Int, q::Float64, pStab::Float64, sat_steps::Int, avg_steps::Int, avg_interval::Int, doPlot::Bool)
    pX = (1. - q)*(1. - pStab)
    pE = q * (1. - pStab)

    # my_stabs = [site_to_op(j,n) for j=1:n]

    my_tableau = initialize_tableau(n)
    # measure_entanglement_scaling(copy(my_tableau), n)
    # for i = 1:2:n
    #     my_tableau = hadamard(my_tableau, i, n)
    #     my_tableau = phase(my_tableau, i, n)
    # end

    measure_entanglement_scaling(copy(my_tableau), n)
    for i = 1:(n-3)
        my_tableau = hadamard(my_tableau, i, n)
        # my_tableau = measure_operator(my_tableau, site_to_op(i,n), n)[1];
        measure_entanglement_scaling(copy(my_tableau), n)
        # print_stabilizer(my_tableau,n)
    end

    for i=1:sat_steps
        my_tableau = update_circuit(my_tableau, n, pStab, pX, pE, i)
        measure_entanglement_scaling(copy(my_tableau), n)
    end

    data_history = zeros(avg_steps, 8)
    # EE, I2, I3, SG, PMX, PMZ, MX, MZ

    SG_matrix = zeros((n ÷ 8), (n ÷ 8))
    PM_X_matrix = zeros((n ÷ 8), (n ÷ 8))
    PM_Z_matrix = zeros((n ÷ 8), (n ÷ 8))

    stab_scaling_history = zeros(avg_steps, n ÷ 2)
    ee_scaling_history = zeros(avg_steps, n)
    mutual_scaling_history = zeros(avg_steps, n)

    step_num = 0
    for i = 1:avg_steps
        data_history[i,1] = compute_bipartite_EE(copy(my_tableau), n)
        data_history[i,2] = measure_antipodal_mutual(copy(my_tableau), n)
        data_history[i,3] = measure_tripartite(copy(my_tableau), n)
        data_history[i,7] = measure_mean_X(copy(my_tableau), n)
        data_history[i,8] = measure_mean_Z(copy(my_tableau), n)
        SG_matrix .+= abs.(measure_spin_glass(copy(my_tableau), n))
        PM_X_matrix .+= abs.(measure_paramagnetic_X(copy(my_tableau), n))
        PM_Z_matrix .+= abs.(measure_paramagnetic_Z(copy(my_tableau), n))
        data_history[i,4] = sum((abs.(SG_matrix ./ i)).^2) / n
        data_history[i,5] = sum((abs.(PM_X_matrix ./ i)).^2) / n
        data_history[i,6] = sum((abs.(PM_Z_matrix ./ i)).^2) / n
        stab_scaling_history[i,:] = measure_stabilizer_scaling(copy(my_tableau), n)
        ee_scaling_history[i,:] = measure_entanglement_scaling(copy(my_tableau), n)
        mutual_scaling_history[i,:] = measure_mutual_scaling(copy(my_tableau), n)
        # (1) Measure scaling of the correlation function for the stabilizer
        # (2) Measure scaling of the entanglement entropy
        for j = 1:avg_interval
            my_tableau = update_circuit(my_tableau, n, pStab, pX, pE, step_num)
            step_num += 1
        end
    end

    SG_matrix ./= avg_steps
    PM_X_matrix ./= avg_steps
    PM_Z_matrix ./= avg_steps

    SG_matrix = (abs.(SG_matrix)).^2
    PM_X_matrix = (abs.(PM_X_matrix)).^2
    PM_Z_matrix = (abs.(PM_Z_matrix)).^2

    SG = sum(SG_matrix) / n
    PM_X = sum(PM_X_matrix) / n
    PM_Z = sum(PM_Z_matrix) / n

    mean_data = mean(data_history,dims=1)[1,:]
    mean_data[4] = SG
    mean_data[5] = PM_X
    mean_data[6] = PM_Z

    mean_ee_scaling = mean(ee_scaling_history, dims=1)[1,:]
    mean_stab_scaling = mean(stab_scaling_history, dims=1)[1,:]
    mean_mutual_scaling = mean(mutual_scaling_history, dims=1)[1,:]

    # if doPlot == true
    #     p = plot(EE_history,lw=1,label="",alpha=.5)
    #     avg_EE_history = cumsum(EE_history) ./ [i for i = 1:avg_steps]
    #     plot!(avg_EE_history, lw=2, label="")
    #     display(p)
    #
    #     p = plot(PM_X_history, lw=2, label="PM X", xscale=:log10)
    #     plot!(p, PM_Z_history, lw=2, label="PM Z", xscale=:log10)
    #     plot!(SG_history, lw=2, label="SG")
    #     display(p)
    #
    # end

    return mean_data, mean_ee_scaling, mean_stab_scaling, mean_mutual_scaling

end

function verify_stabilizer(tableau, n)
    for i = (n+1):(2n)
        for j = (n+1):(2n)
            if check_commute(tableau[i,:], tableau[j,:], n) == false
                return false
            end
        end
    end
    return true
end

function print_stabilizer(tableau, n)
    my_stabs = tableau[(n+1):(2n), 1:(2n)]
    strs = ["" for i=1:n]
    for i=1:n
        for j = 1:n
            if my_stabs[i,j] == 1 && my_stabs[i,j+n] == 1
                strs[i] = string(strs[i], "Y")
            elseif my_stabs[i,j] == 1 && my_stabs[i,j+n] == 0
                strs[i] = string(strs[i], "X")
            elseif my_stabs[i,j+n] == 1
                strs[i] = string(strs[i],"Z")
            else
                strs[i] = string(strs[i],"I")
            end
        end
        println(strs[i])
    end
end

function main(;q=0,n=16, num_reps=2, avg_steps = 100)
    # p_range = [i for i=0:0.1:1]
    # p_range = [0.4, 0.45, 0.5, 0.55, 0.6]
    p_range = [0.4, 0.5, 0.6, 0.7]
    sat_steps = 200
    avg_interval = 31

    all_dat = zeros(length(p_range), num_reps, 8)
    all_ee_scaling = zeros(length(p_range), num_reps, n)
    all_stab_scaling = zeros(length(p_range), num_reps, n ÷ 2)
    all_mutual_scaling = zeros(length(p_range), num_reps, n)
    # Order is EE, I2, I3, SG, PMX, PMZ, MX, MZ

    for i = 1:length(p_range)
        println("\n", i)
        for j = 1:num_reps
            if (j % 10 == 0)
                print(".")
            end
            res_mean, res_ee, res_stab, res_mut = run_circuit(n,q,p_range[i], sat_steps, avg_steps, avg_interval, false)
            all_dat[i,j,:] .= res_mean
            all_ee_scaling[i,j,:] .= res_ee
            all_stab_scaling[i,j,:] .= res_stab
            all_mutual_scaling[i,j,:] .= res_mut
            # all_dat[i,j,:] .= run_circuit(n,q,p_range[i], sat_steps, avg_steps, avg_interval, false)
        end
    end
    println("\n")

    fname = string("513_code/code513_q", round(q,digits=2), "_N", n, "_data_pure_3.txt")
    mean_dat = mean(all_dat, dims=2)[:,1,:]
    std_dat = std(all_dat, dims=2)[:,1,:]
    to_write = cat(p_range, mean_dat, std_dat, dims=2)
    open(fname, "w") do io
        writedlm(io, to_write, ',')
    end
    fname_ee = string("513_code/code513_q", round(q,digits=2), "_N", n, "_EE_data_pure_3.txt")
    fname_stab = string("513_code/code513_q", round(q,digits=2), "_N", n, "_stab_data_pure_3.txt")
    fname_mut = string("513_code/code513_q", round(q,digits=2), "_N", n, "_I2_data_pure_3.txt")
    mean_ee = mean(all_ee_scaling, dims=2)[:,1,:]
    std_ee = std(all_ee_scaling, dims=2)[:,1,:]
    to_write_ee = cat(mean_ee, std_ee, dims=2)
    mean_stab = mean(all_stab_scaling, dims=2)[:,1,:]
    std_stab = std(all_stab_scaling, dims=2)[:,1,:]
    to_write_stab = cat(mean_stab, std_stab, dims=2)
    mean_mut = mean(all_mutual_scaling, dims=2)[:,1,:]
    std_mut = std(all_mutual_scaling, dims=2)[:,1,:]
    to_write_mut = cat(mean_mut, std_mut, dims=2)
    open(fname_ee, "w") do io_ee
        writedlm(io_ee, to_write_ee, ',')
    end
    open(fname_stab, "w") do io_stab
        writedlm(io_stab, to_write_stab, ',')
    end
    open(fname_mut, "w") do io_mut
        writedlm(io_mut, to_write_mut, ',')
    end

    p = plot(p_range, mean_dat[:,4] .* (64/n), yerror = std_dat[:,4] .* (64/n), markershape=:circle, label="SG")
    plot!(p_range, mean_dat[:,5] .* (64/n), yerror = std_dat[:,5] .* (64/n), markershape=:circle, label="PM X")
    plot!(p_range, mean_dat[:,6] .* (64/n), yerror = std_dat[:,6] .* (64/n), markershape=:circle, label="PM Z")
    plot!(p_range, mean_dat[:,7], yerror = std_dat[:,7], markershape=:circle, label="MX")
    plot!(p_range, mean_dat[:,8], yerror = std_dat[:,8], markershape=:circle, label="MZ")
    Plots.savefig(p,string("513_code/code513_q", round(q,digits=2), "_N", n, "_order_param_pure_3.png"))

    p2 = plot(p_range, -mean_dat[:,1], yerror=std_dat[:,1], markershape=:circle, label="")
    Plots.savefig(p2,string("513_code/code513_q", round(q,digits=2), "_N", n, "_EE_pure_3.png"))

    return 1
end

# n = 16
# if length(ARGS) >= 1
#     n = parse(Int64, ARGS[1])
# end
# @show n
#
# for q=0.:.1:1
#     println("q=$q")
#     main(q=q, n=n, num_reps=10)
#     break
# end

println("loaded")
