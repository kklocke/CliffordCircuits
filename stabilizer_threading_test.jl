using Distributed

global n = 18
if length(ARGS) >= 1
    n = parse(Int64, ARGS[1])
end
@show n

# @distributed for i in 1:2
#     r = test_func(i) # rand()
#     r2 = site_to_op_test(5)
#     println(i," ",r," ",size(r2)[1])
#     # main(q=0.1*(i-1), n=16, avg_steps=100, num_reps=3)
# end

println(size(site_to_op_test(4))[1])

# @distributed for i in 1:2
#     r = size(site_to_op_test(5))
#     println(i, " ", r)
#     main(q=0.1*(i-1), n=16, avg_steps=100, num_reps=3)
# end

# addprocs(3);
println("num procs = ", nprocs())

xList = [i for i in 1:11]
# xMean = pmap(site_to_op_test, xList)
# xMean = pmap(x->site_to_op_test(x), xList)
xMean = pmap(x->main(q=0.1*(x-1), n=n, num_reps=2, avg_steps=150), xList);
println(xMean)
