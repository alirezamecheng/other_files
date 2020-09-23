# @info("Loading Plots Package")
# using Plots
@info("Loading FDM 2D module")
include("modules/FDM_D2Q4.jl")
using .FDM_D2Q4
@info("Loading LBM BGK D2Q4 module")
include("modules/BGK_D2Q4.jl")
using .LBM_BGK_D2Q4
include("modules/BGK_D2Q4_Parallel.jl")
using .LBM_BGK_D2Q4_Parallel
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

    # timeMax = [.1 .15 .2]
    # tauSet = [1 2 3 4 5 6 7 8 9 10]
    timeMax = [.2]
    tauSet = [2]
    ### timing output ###
    if isfile("benchmarking.csv")
        timingIO = open("benchmarking.csv","a") # ==> فایل خروجی نتایج زمان بندی
    else
        timingIO = open("benchmarking.csv","a") # ==> فایل خروجی نتایج زمان بندی
        println(timingIO,"threadNumber,tstar,tau,LBM_evals,FDM_evals,LBMTime_median,LBMBenchParallel_median,LBMTime_mean,LBMBenchParallel_mean") # ==> فاکتور هایی که در فایل زمان بندی ذخیره می شوند. 
    end

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
            LBMBenchParallelSteps = 4.0 * tstar * (N_f-1)^2 # number of time steps for FDM

            println("Time = ", tstar)
            println("LBMTimeSteps = ",LBMTimeSteps)
            println("LBMBenchParallelSteps = ", LBMBenchParallelSteps)
            # @info("Press enter to procede: ")
            # think = readline()

            # @info("Generating the FDM Problem")
            # FDM_Problem = FDM_D2Q4.FDM_2D_Problem_maker(L = L, alpha= alpha, N = N_f, 
            # max_time = Int64(round(FDMTimeSteps)), Tinitial= Tinitial, T_top = T_top, T_left = T_left,
            # T_down = T_down, T_right = T_right)

            @info("Generating the LBM Problem parallel")
            LBM_Problem_Parallel = LBM_BGK_D2Q4_Parallel.LBM_BGK_D2Q4_Problem_maker(N = N_l, tau = tau, 
            max_time = Int64(round(LBMTimeSteps)), Tinitial= Tinitial, T_top = T_top, T_left = T_left, 
            T_down = T_down, T_right = T_right)

            @info("Generating the LBM Problem")
            LBM_Problem = LBM_BGK_D2Q4.LBM_BGK_D2Q4_Problem_maker(N = N_l, tau = tau, 
            max_time = Int64(round(LBMTimeSteps)), Tinitial= Tinitial, T_top = T_top, T_left = T_left, 
            T_down = T_down, T_right = T_right)

            @info("Solving the Problems")
            LBMbench = @benchmarkable T_LBM = LBM_BGK_D2Q4.solve($LBM_Problem)
            LBMBenchParallel = @benchmarkable T_LBM = LBM_BGK_D2Q4_Parallel.solve($LBM_Problem_Parallel)
            # FDMbench = @benchmarkable T_FDM = FDM_D2Q4.solve($FDM_Problem)
            LBMTime = run(LBMbench,seconds=10)
            LBMBenchParallel = run(LBMBenchParallel,seconds=10)
            println(timingIO,Threads.nthreads(),",",tstar,",",tau,",",size(LBMTime.times,1),",",size(LBMBenchParallel.times,1),",",median(LBMTime.times),",",median(LBMBenchParallel.times),",",mean(LBMTime.times),",",mean(LBMBenchParallel.times))
            flush(timingIO)
        end #for tau
    end # end for timeMax