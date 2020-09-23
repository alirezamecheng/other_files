# @info("Loading Plots Package")
# using Plots
@info("Loading FDM 2D module")
include("modules/FDM_D2Q4.jl")
using .FDM_D2Q4
@info("Loading LBM BGK D2Q4 module")
include("modules/BGK_D2Q4.jl")
using .LBM_BGK_D2Q4
@info("Loading Benchmarking tooks")
using BenchmarkTools
@info("Loading Statistics package")
using Statistics
@info("Loading tools")
include("tools/tools.jl")


# function julia_main()
    
    ii = 1

    # Defining the Problem
    L = 1.0 # length [m]
    N = 201 # number of nodes
    # tstar = 1 # dimention less solution time
    # tau = 2.1 # relaxation time for LBM # max 5.1 is good
    alpha = 5.0*10^-6 # diffusivity for FDM [For mercury]

    # Dricklet boundary condition Temperature in [C]
    T_top = 1
    T_left = 1
    T_down = 1
    T_right = 1

    # initial condition 
    Tinitial = 0

    N_f = N # number of nodes for FDM
    N_l = N # number of nodes for LBM

    dx = L/(N-1) # FDM dx
    x = 0:dx:L     # x domain
    ### Allocation of single output 
    T_FDM = Array{Float64,2}(undef,N_f,N_f)
    T_LBM = Array{Float64,2}(undef,N_l,N_l)

    # timeMax = [.01 .05 .1 .15 .2 .25 .3 .35 .4]
    # tauSet = [1 2 3 4 5 6 7 8 9 10]
    timeMax = [.1 .15 .2]
    tauSet = [1 2 3 4 5 6 7 8 9 10]
    ### timing output ###
    timingIO = open("benchmarking.csv","w")
    println(timingIO,"tstar,tau,LBM_evals,FDM_evals,LBMTime_median,FDMTime_median,LBMTime_mean,FDMTime_mean")
    for tstar in timeMax

        for tau in tauSet
            println("========= tau = $(tau) =========")
            println("======== time = $(tstar) ========")
        
            dx_f = L/(N_f-1) # dx for FDM
            dt_f = (1/4)*(dx_f^2)/(alpha) # dt for FDM based on [Fo <= 1/4] Which has the maximum time step

            dt_l = 1.0 # LBM time step in lattice units 
            dx_l = 1.0 # LBM  space step in lattice units

            # tau = 2*alpha*dt_l / dx_l^2 + 0.5
            alpha_l = ((dx_l)^2/(2*dt_l))*(tau - 0.5) # diffusivity for LBM in Lattice units
            # LBMTimeSteps = tstar * N_l^2 * dx_l^2 / (alpha_l * dt_l)
            LBMTimeSteps = 2.0 * tstar * (N_l-1)^2 /(tau - 0.5) # number of time steps for LBM
            FDMTimeSteps = 4.0 * tstar * (N_f-1)^2 # number of time steps for FDM

            println("Time = ", tstar)
            println("LBMTimeSteps = ",LBMTimeSteps)
            println("FDMTimeSteps = ", FDMTimeSteps)
            # @info("Press enter to procede: ")
            # think = readline()

            @info("Generating the FDM Problem")
            FDM_Problem = FDM_D2Q4.FDM_2D_Problem_maker(L = L, alpha= alpha, N = N_f, 
            max_time = Int64(round(FDMTimeSteps)), Tinitial= Tinitial, T_top = T_top, T_left = T_left,
            T_down = T_down, T_right = T_right)

            @info("Generating the LBM Problem")
            LBM_Problem = LBM_BGK_D2Q4.LBM_BGK_D2Q4_Problem_maker(N = N_l, tau = tau, 
            max_time = Int64(round(LBMTimeSteps)), Tinitial= Tinitial, T_top = T_top, T_left = T_left, 
            T_down = T_down, T_right = T_right)

            @info("Solving the Problems")
            LBMbench = @benchmarkable T_LBM = LBM_BGK_D2Q4.solve($LBM_Problem)
            FDMbench = @benchmarkable T_FDM = FDM_D2Q4.solve($FDM_Problem)
            LBMTime = run(LBMbench,seconds=10)
            FDMTime = run(FDMbench,seconds=10)
            println(timingIO,tstar,",",tau,",",size(LBMTime.times,1),",",size(FDMTime.times,1),",",median(LBMTime.times),",",median(FDMTime.times),",",mean(LBMTime.times),",",mean(FDMTime.times))
            flush(timingIO)
        end #for tau
    end # end for timeMax

    # @info("Saving the outputs to files")
    # for i = 1:size(timeMax,2)
    #     saveTP(T_FDM[i],"FDM$(timeMax[i])")
    #     saveTP(T_LBM[i],"LBM$(timeMax[i])")
    # end
    # for i = 1:size(timeMax,2)
    #     saveCenterLines(T_FDM[i],"FDM$(timeMax[i])")
    #     saveCenterLines(T_LBM[i],"LBM$(timeMax[i])")
    # end
    # @info("You can find the tp files in the results folder")

    # @info("Compiling Plots Package")
    # # pyplot()
    # i = 1
    # p = plot(x,T_FDM[i][Int64(floor(N_f/2)),:], label = "FDM t = $(timeMax[i])")
    # # scatter!(p, x[1:10:end],T_LBM[i][Int64(floor(N_l/2)),1:10:end], label = "LBM t = $(timeMax[i])")
    # for i = 2:size(timeMax,2)
    #     plot!(p ,x,T_FDM[i][Int64(floor(N_f/2)),:], label = "FDM t = $(timeMax[i])")
    #     # scatter!(p, x[1:10:end],T_LBM[i][Int64(floor(N_l/2)),1:10:end], label = "LBM t = $(timeMax[i])")
    # end

    # @info("Plotting the results")
    # yaxis!("Temperature [°C]") 
    # xaxis!("Length [m]") 
    # plot!(p,legendtitle="Solution time [s]")
    # plot!(p,legend=false)
    # savefig(p,"RESULTS.pdf")
    
    # @info("Results are saved.")
    # if Sys.islinux()
    #     run(Cmd(`xdg-open RESULTS.pdf`))
    # end

    # return T_LBM, T_FDM
# end # end function julia_main()
# T_LBM, T_FDM = julia_main();

