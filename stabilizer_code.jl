include("stabilizer_helper.jl")

function measure_spin_glass(tableau, n)
    A = [i for i=1:(n÷8)] .+ (n÷4)
    B = A .+ (n÷2)
    my_matrix = zeros(length(A), length(B))
    for i in A
        for j in B
            t,z12,ez = measure_zz(copy(tableau), i, j, n)
            my_matrix[i - A[1] + 1,j - B[1] + 1] = ez
        end
    end
    return my_matrix
end

function measure_paramagnetic(tableau, n)
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

function update_circuit(tableau, n, pZZ, pX, pE, i)
    if i % 2 == 0
        for j = 1:n
            next_site = j+1
            if j == n
                next_site = 1
            end
            if rand() < pZZ
                tableau, z12, ez = measure_zz(tableau, j, next_site, n)
            end
        end
    else
        for j = 1:n
            r = rand()
            if r < pX
                tableau, x1, ex = measure_x(tableau, j, n)
            elseif r < (pE + pX)
                tableau = simple_dephase(tableau, j, n)
            end
        end
    end
    return tableau
end

function run_circuit(n, q, pZZ, sat_steps, avg_steps, avg_interval, doPlot)
    pX = (1. - q)*(1. - pZZ)
    pD = q * (1. - pZZ);

    my_tableau = initialize_tableau(n);
    for i = 1:n
        my_tableau = hadamard(my_tableau, i, n)
    end

    for i = 1:sat_steps
        my_tableau = update_circuit(my_tableau, n, pZZ, pX, pD, i)
    end

    EE_history = zeros(avg_steps)
    SG_history = zeros(avg_steps)
    PM_history = zeros(avg_steps)
    PM_matrix = zeros((n ÷ 8), (n ÷ 8))
    SG_matrix = zeros((n ÷ 8), (n ÷ 8))

    step_num = 0
    for i = 1:avg_steps
        for j = 1:avg_interval
            my_tableau = update_circuit(my_tableau, n, pZZ, pX, pD, step_num)
            step_num += 1
        end
        EE_history[i] = compute_bipartite_EE(my_tableau, n)
        SG_matrix .+= abs.(measure_spin_glass(my_tableau, n))
        PM_matrix .+= abs.(measure_paramagnetic(my_tableau, n))
        SG_history[i] = sum((abs.(SG_matrix ./ i)).^2) / n
        PM_history[i] = sum((abs.(PM_matrix ./ i)).^2) / n
    end

    SG_matrix ./= avg_steps
    PM_matrix ./= avg_steps

    SG_matrix = (abs.(SG_matrix)).^2
    PM_matrix = (abs.(PM_matrix)).^2

    SG = sum(SG_matrix) / n
    PM = sum(PM_matrix) / n

    EE = mean(EE_history)

    if doPlot == true
        p = plot(EE_history,lw=1,label="",alpha=.5)
        avg_EE_history = cumsum(EE_history) ./ [i for i = 1:avg_steps]
        plot!(avg_EE_history, lw=2, label="")
        display(p)

        p = plot(PM_history, lw=2, label="PM", xscale=:log10)
        plot!(SG_history, lw=2, label="SG")
        display(p)
    end

    return (EE, SG, PM)
end

function main()
    p_range = [i for i = 0:0.1:1]
    sat_times = [1000,1000,1000,1000,1000,1500,1500,1500,5000,5000,5000]
    n = 32
    q = 0.5
    sat_steps = 300
    avg_steps = 1000
    avg_interval = 31

    num_reps = 2

    sg_vals = zeros(length(p_range), num_reps)
    pm_vals = zeros(length(p_range), num_reps)
    ee_vals = zeros(length(p_range), num_reps)

    for i = 1:length(p_range)
        println("\n", i)
        for j = 1:num_reps
            if (j % 10 == 0)
                print(".")
            end
            ee, sg, pm = run_circuit(n,q,p_range[i], sat_steps, avg_steps, avg_interval, false)
            ee_vals[i,j] += ee
            sg_vals[i,j] += sg
            pm_vals[i,j] += pm
        end
    end
    println("\n")

    p = plot(p_range, mean(pm_vals,dims=2)  .* (64 / n), yerror=std(pm_vals, dims=2) .* (64/n), markershape=:circle,label="PM")
    plot!(p_range, mean(sg_vals,dims=2) .* (64 / n), yerror=std(sg_vals,dims=2) .* (64/n), markershape=:circle,label="SG",legend=:outerright)
    Plots.savefig(p,string("li_fisher_qPt5_N32_order_param.png"))

    p2 = plot(p_range, -mean(ee_vals,dims=2), yerror=std(ee_vals,dims=2), markershape=:circle, label="")
    Plots.savefig(p2,string("li_fisher_qPt5_N32_EE.png"))
end

main()
