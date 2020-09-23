# version 3 
#= 
Equations from: 
Kałuża, G. (2012). 
The Numerical Solution of Richards Equation Using the Lattice Boltzmann Method. 
Applied Mechanics and Materials, 11(1), 23–30. Retrieved 
from http://www.amcm.pcz.pl/

key \cite{Kauza2012}
=#


module D1Q2

const Q = 2
const W = [1 / 2, 1 / 2]


struct  LBM_D1Q2_Heat_conduction_Problem
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
    f
end 

function initialCondition!(f::Array{Float64,2}, Tinitial::Float64)
    for k = 1:Q
        @simd for i = 1:size(f, 1)
            @inbounds f[i,k] = W[k] * Tinitial
        end
    end
    return f
end

function collision!(f::Array{Float64,2}, 
                    fpost::Array{Float64,2}, 
                    T::Array{Float64,1}, 
                    tau::Float64, 
                    dt::Float64)
    for k = 1:Q
        @simd for i = 1:size(fpost, 1)
            @inbounds fpost[i,k] = ((1 - (dt / tau)) * f[i,k]) + (dt / tau) * W[k] * T[i]
        end
    end
    return fpost
end


function streaming!(fpost::Array{Float64,2})
    @simd for j = (size(fpost, 1) - 1):-1:1 # should be this way otherwise the method will not work
        @inbounds fpost[j + 1,1] = fpost[j,1]
    end
    @simd for j = 2:size(fpost, 1)
        @inbounds fpost[j - 1,2] = fpost[j,2]
    end
    return fpost
end
function boundaryCondition!(fpost::Array{Float64,2},
                            Tw::Float64,
                            Te::Float64;
                            landa::Float64,
                            alpha::Float64,
                            h::Float64)
    @inbounds fpost[1,1] = Tw - fpost[1,2]
    @inbounds fpost[end,2] = (landa/(landa + alpha*h))*(fpost[end-1,1] + fpost[end-1,2]) - fpost[end,1] + (alpha*h/(landa+alpha*h))*Te
    
    return fpost
end

function calcT!(f::Array{Float64,2},T::Array{Float64,1})
    @simd for i = 1:size(f, 1)
        @inbounds T[i] = f[i,1] + f[i,2]
    end
    return T
end
function timeMarching!(f::Array{Float64,2}, fpost::Array{Float64,2})
    for j = 1:Q
        @simd for i = 1:size(fpost,1)
            @inbounds f[i,j] = fpost[i,j]
        end
    end
    return f
end

function LBM_Heat_Problem_Generator(
                                    ;
                                    L::Real,
                                    h::Real,
                                    max_time::Real,
                                    dt::Real,
                                    landa::Real,
                                    c::Real,
                                    Tw::Real,
                                    Te::Real,
                                    Tinitial::Real)

    return   LBM_D1Q2_Heat_conduction_Problem(
                                                Float64(L),
                                                Float64(h), 
                                                Float64(max_time), 
                                                Float64(dt), 
                                                Float64(landa), 
                                                Float64(c), 
                                                Float64(Tw), 
                                                Float64(Te), 
                                                Float64(Tinitial))
end # end of [LBM_Heat_Problem_Generator] 



function solver(P::LBM_D1Q2_Heat_conduction_Problem)
    Nx = Int64((P.L) / (P.h) + 1)
    max_time = P.max_time
    c = P.c
    alpha = (P.landa) / (c)
    dt = P.dt
    nu = (P.h) / (P.dt)
    Tw = P.Tw
    Te = P.Te
    Tinitial = P.Tinitial

    # Allocation 

    f = Array{Float64,2}(undef, Nx, Q)
    fpost = Array{Float64,2}(undef, Nx, Q)
    T = Array{Float64,1}(undef, Nx)

    # special variables for methods
    tau = (alpha / (nu^2)) + (dt / 2)

    f = initialCondition!(f, Tinitial)
    T .= Tinitial
    
    for iteration = 0:dt:max_time
        # println("Iteration : ", iteration)
        fpost = collision!(f, fpost, T, tau, dt)
        fpost = streaming!(fpost)
        fpost = boundaryCondition!(fpost, Tw, Te,landa = P.landa, alpha = alpha, h =P.h)
        T = calcT!(fpost,T)
        f = timeMarching!(f, fpost)
    end


    return solution(T, f)
end 

export LBM_Heat_Problem_Generator, solver
end # End of [D1Q2] module 

# using .D1Q2 
# 
# function julia_main()
#     problem = D1Q2.LBM_Heat_Problem_Generator(
#                                         L = 0.05,
#                                         h = 0.00125, 
#                                         max_time = 120,
#                                         dt = 0.1,
#                                         landa = 35,
#                                         c = 4.857*10^6,
#                                         Tw = 1,
#                                         Te = 1,
#                                         Tinitial = 0.1)
#     solution = D1Q2.solver(problem)
# end
# 
# @time julia_main()
