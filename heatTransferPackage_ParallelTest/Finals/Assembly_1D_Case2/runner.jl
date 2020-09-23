include("./tools/packagesAndModulesLoader.jl") # including Tools and Packages

function runs()
    L = 0.05
    h = 0.00125
    max_time = 120
    dt = 0.001
    landa = 35
    c = 4.875 * 10^6
    Tw = 150
    Te = 20
    Tinitial = 0

    dx = h
    alpha = landa / c
    Nx = Int64(L / dx + 1)

    # checking the divergence criteria for FDM
    if (alpha*dt) > (1/2*dx)
        @warn("The convergence criteria is not satisfied the problem will diverge.")
        println("inputs are: ")
        println("L = ",L)
        println("dx = ",dx)
        println("dt = ",dt)
        println("α = ",alpha) 
        println("Nx = ",Nx)
        println("Fo = ",alpha*dt/dx)
    else 
        @info("Solving with the following inputs:")
        println("L = ",L)
        println("dx = ",dx)
        println("dt = ",dt)
        println("Max_TIME = ",max_time)
        println("α = ",alpha) 
        println("Nx = ",Nx)
        println("Fo = ",alpha*dt/dx)
    end

    FDM_Solution = []
    LBM_Solution = []
    ii = 1
    solutionTimes = [2 6 10 40 120]
    for max_time in solutionTimes
        @info("Solving for maxTime  = $(max_time)")
        problem1 = FDMHeatTransfer.FDM_1D_convection_Heat_Transfer_Problem_Generator(
                                                                                    L = L,
                                                                                    h = h,
                                                                                    max_time = max_time,
                                                                                    dt = dt,
                                                                                    landa = landa,
                                                                                    c = c,
                                                                                    Tw = Tw,
                                                                                    Te = Te,
                                                                                    Tinitial = Tinitial)
        @time  push!(FDM_Solution,FDMHeatTransfer.solve(problem1))

        problem2 = D1Q2.LBM_Heat_Problem_Generator(
                                                    L = L,
                                                    h = h,
                                                    max_time = max_time,
                                                    dt = dt,
                                                    landa = landa,
                                                    c = c,
                                                    Tw = Tw,
                                                    Te = Te,
                                                    Tinitial = Tinitial)
        @time push!(LBM_Solution,D1Q2.solver(problem2))
        ii = ii + 1
        @info("Solving for maxTime  = $(max_time) is finished.")
    end


    x = 0:dx:L



    @info("Ploting the results")
    i = 1
    # different style plot
    # fontsize = 10
    # p = plot(x,FDM_Solution[1].T,label = "FDM t = $(solutionTimes[1])",legend=:topright,legendfontsize=6,yguidefontsize=fontsize,xtickfontsize=fontsize,ytickfontsize=fontsize,xguidefontsize=fontsize,titlefontsize=fontsize)
    # # plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])")
    # scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])",markershape=:rect)
    # i = i + 1
    # plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])",linestyles = :dash)
    # scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])",markershape=:pentagon)
    # i = i + 1
    # plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])",linestyles = :dot)
    # scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])",markershape=:hexagon)
    # i = i + 1
    # plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])",linestyles = :dashdot)
    # scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])",markershape=:heptagon)
    # i = i + 1
    # plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])",linestyles = :dashdotdot)
    # scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])",markershape=:octagon)
    # # title!("Temperature distribution in the wall")
    # yaxis!("Temperature [°C]") 
    # xaxis!("Length [m]") 
    # plot!(p,legendtitle="Solution times [s]")

    fontsize = 10
    p = plot(x,FDM_Solution[1].T,label = "FDM t = $(solutionTimes[1])",legend = :topright,legendfontsize=6,yguidefontsize=fontsize,xtickfontsize=fontsize,ytickfontsize=fontsize,xguidefontsize=fontsize,titlefontsize=fontsize)
    # plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])")
    scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])")
    i = i + 1
    plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])")
    scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])")
    i = i + 1
    plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])")
    scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])")
    i = i + 1
    plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])")
    scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])")
    i = i + 1
    plot!(p,x,FDM_Solution[i].T,label = "FDM t = $(solutionTimes[i])")
    scatter!(x,LBM_Solution[i].T,label = "LBM t = $(solutionTimes[i])")
    # title!("Temperature distribution in the wall")
    yaxis!("Temperature [°C]") 
    xaxis!("Length [m]") 
    plot!(p,legendtitle="Solution times [s]")




    savefig(p,"RESULTS.pdf")
    @info("Results are saved.")
    if Sys.islinux()
        run(Cmd(`xdg-open RESULTS.pdf`))
    end
end

function timings()
    L = 0.05
    h = 0.00125
    max_time = 120
    dt = 0.001
    landa = 35
    c = 4.875 * 10^6
    Tw = 0
    Te = 100
    Tinitial = 0

    dx = h
    alpha = landa / c
    Nx = Int64(L / dx + 1)
    @info("Time Benchmarking The Methods")
    max_time = 10
    problem10 = FDMHeatTransfer.FDM_1D_convection_Heat_Transfer_Problem_Generator(L = L,
                                                                                h = h,
                                                                                max_time = max_time,
                                                                                dt = dt,
                                                                                landa = landa,
                                                                                c = c,
                                                                                Tw = Tw,
                                                                                Te = Te,
                                                                                Tinitial = Tinitial)
    FDMTimes = @benchmark FDMHeatTransfer.solve($(problem10))
    # FDMHeatTransfer.solve(problem10)

    problem20 = D1Q2.LBM_Heat_Problem_Generator(L = L,
                                                h = h,
                                                max_time = max_time,
                                                dt = dt,
                                                landa = landa,
                                                c = c,
                                                Tw = Tw,
                                                Te = Te,
                                                Tinitial = Tinitial)
    LBMTimes = @benchmark D1Q2.solver($(problem20))

    @info("Preparing benchmarking time result")
    println("Minimum time for FDM in ",size(FDMTimes.times,1)," repeats is : ", minimum(FDMTimes.times)) 
    println("Minimum time for LBM in ",size(LBMTimes.times,1)," repeats is : ", minimum(LBMTimes.times))
end





function julia_main()
    runs()
    println("\n\n")
    # timings()
end
julia_main()
