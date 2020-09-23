# finite difference method for 1d heat conduction
#= from [Kauza2012]
    explicit scheme of the finite difference method (FDM) [6-8]
=#

module FDMHeatTransfer

struct FDM_1D_convection_Heat_Transfer_Problem
    L::Float64
    h::Float64
    max_time::Float64
    dt::Float64
    landa::Float64
    c::Float64
    Tw::Float64
    Te::Float64
    Tinitial::Float64
end
mutable struct solution
    T
end
function stencile!(;Tnew::Array{Float64,1}, T::Array{Float64,1}, alpha::Float64, dt::Float64, h::Float64)
    for i = 2:(size(T, 1)-1)
        @inbounds Tnew[i] = (1 - (2 * alpha * dt) / (h^2)) * T[i] + (alpha * dt) / (h^2) * (T[i - 1] + T[i + 1])
    end
    return Tnew
end
function initialCondition!(;T::Array{Float64,1}, Tinitial::Float64)
    for i = 1:size(T, 1)
        @inbounds T[i] = Tinitial
    end
    return T
end
function stepForward!(;Tnew::Array{Float64,1}, T::Array{Float64,1})
    for i = 1:size(T, 1)
        @inbounds T[i] = Tnew[i]
    end
    return T
end
function boundaryCondition!(;Tnew::Array{Float64,1}, 
    T::Array{Float64,1}, Tw::Float64, Te::Float64,
    alpha::Float64, h::Float64, landa::Float64, Ta::Float64)

    @inbounds Tnew[1] = Tw 
    # @inbounds T[end] = Te
    @inbounds Tnew[end] = landa/(landa + alpha*h) *Tnew[end-1] + alpha*h/(landa + alpha*h)*Ta 
    return T 
end
function FDM_1D_convection_Heat_Transfer_Problem_Generator(;
                                                            L::Real,
                                                            h::Real,
                                                            max_time::Real,
                                                            dt::Real,
                                                            landa::Real,
                                                            c::Real,
                                                            Tw::Real,
                                                            Te::Real,
                                                            Tinitial::Real)
    return FDM_1D_convection_Heat_Transfer_Problem(Float64(L),
                                                    Float64(h),
                                                    Float64(max_time),
                                                    Float64(dt),
                                                    Float64(landa),
                                                    Float64(c),
                                                    Float64(Tw),
                                                    Float64(Te),
                                                    Float64(Tinitial))
end
function solve(P::FDM_1D_convection_Heat_Transfer_Problem)
    L = P.L
    dx = P.h
    MAX_TIME = P.max_time
    dt = P.dt
    landa = P.landa
    c = P.c
    Tw = P.Tw
    Te = P.Te
    Tinitial = P.Tinitial

    # Ta = 20.0
    alpha = landa / c
    Nx = Int64(L / dx + 1)

    # if (alpha*dt) > (1/2*dx)
    #     @warn("The convergence criteria is not satisfied the problem will diverge.")
    #     println("inputs are: ")
    #     println("L = ",L)
    #     println("dx = ",dx)
    #     println("dt = ",dt)
    #     println("α = ",alpha) 
    #     println("Nx = ",Nx)
    #     println("Fo = ",alpha*dt/dx)
    # else 
    #     @info("Solving with the following inputs:")
    #     println("L = ",L)
    #     println("dx = ",dx)
    #     println("dt = ",dt)
    #     println("Max_TIME = ",MAX_TIME)
    #     println("α = ",alpha) 
    #     println("Nx = ",Nx)
    #     println("Fo = ",alpha*dt/dx)
    # end

    # Allocations
    T = Array{Float64,1}(undef, Nx)
    Tnew = Array{Float64,1}(undef, Nx)

    T = initialCondition!(T = T, Tinitial = Tinitial)
    T = boundaryCondition!(Tnew = Tnew, T = T, Tw = Tw, Te = Te, alpha = alpha, landa = landa, h = dx, Ta = Te)

    for i = 0:dt:MAX_TIME
        Tnew = stencile!(Tnew = Tnew, T = T, alpha = alpha, dt = dt, h = dx)
        Tnew = boundaryCondition!(Tnew = Tnew, T = T, Tw = Tw, Te = Te, alpha = alpha, landa = landa, h = dx, Ta = Te)
        T = stepForward!(Tnew = Tnew, T = T)
    end

    return solution(T)
end
export solve, FDM_1D_convection_Heat_Transfer_Problem_Generator
end # end module [FDMHeatTransfer]

# using .FDMHeatTransfer
# 
# 
# function julia_main()
#     problem = FDMHeatTransfer.FDM_1D_convection_Heat_Transfer_Problem_Generator(
#     L=0.05,
#     h=0.00125,
#     max_time=120,
#     dt=0.1,
#     landa=35,
#     c=4.875*10^6,
#     Tw=0,
#     Te=100,
#     Tinitial=0)
#     solution = FDMHeatTransfer.solve(problem)
# end
# julia_main()
