# math from [mohammad]
module FDM_D2Q4
export solve, FDM_2D_Problem_maker
struct FDM_2D_Problem
    L::Float64
    alpha::Float64
    N::Int64
    max_time::Int64
    Tinitial::Float64
    T_top::Float64
    T_left::Float64
    T_down::Float64
    T_right::Float64
end

function FDM_2D_Problem_maker(;L::Number,
                                alpha::Number,
                                N::Int64,
                                max_time::Int64,
                                Tinitial::Number,
                                T_top::Number,
                                T_left::Number,
                                T_down::Number,
                                T_right::Number)
    P = FDM_2D_Problem(Float64(L),
                        Float64(alpha),
                        Int64(N),
                        Int64(max_time),
                        Float64(Tinitial),
                        Float64(T_top),
                        Float64(T_left),
                        Float64(T_down),
                        Float64(T_right)
                        )

    return P
end

function body!(T::Array{Float64}, Tnew::Array{Float64};alpha::Float64, dt::Float64, dx::Float64, dy::Float64)
    for j = 2:size(T,2)-1
        for i = 2:size(T,1)-1
            termX = (T[i+1,j] + T[i-1,j])/(dx^2)
            termY = (T[i,j+1] + T[i,j-1])/(dy^2)
            dd = 1/(dx^2) + 1/(dy^2)
            Tnew[i,j] = T[i,j] + dt*alpha * (termX + termY - 2.0*T[i,j] * dd)
        end
    end
    return Tnew
end
function BC!(T::Array{Float64}; T_top::Float64, T_down::Float64, T_left::Float64, T_right::Float64)
    Nx = size(T,1)
    Ny = size(T,2)
    # T_top
    j = Ny
    for i = 1:Nx
        T[i,j] = T_top
    end
    # T_down
    j = 1
    for i = 1:Nx
        T[i,j] = T[i,j+1]
    end
    # T_left
    i = 1
    for j = 1:Ny
        T[i,j] = T_left
    end
    # T_right
    i = Nx
    for j = 1:Ny
        T[i,j] = T_right
    end
    return T  
end
function timeMarche!(T::Array{Float64}, Tnew::Array{Float64})
    for j = 1:size(T,2)
        for i = 1:size(T,1)
            T[i,j] = Tnew[i,j]
        end
    end
    return T
end
function initialize!(T::Array{Float64}, Tinitail::Float64)
    return T .=Tinitail
end
function solve(P::FDM_2D_Problem)
    # parameters
    Tinitial = P.Tinitial
    T_top = P.T_top
    T_left = P.T_left
    T_down = P.T_down
    T_right = P.T_right

    L = P.L
    alpha = P.alpha
    # calc other parameters
    N = P.N
    Nx = N
    Ny = N 

    dx = L/(N-1)
    dy = dx
    dt = (1/4) * (L^2 /((N-1)^2 * alpha))
    max_time = P.max_time

    # Allocation
    T = Array{Float64,2}(undef, Nx, Ny)
    Tnew = Array{Float64,2}(undef, Nx, Ny)

    # initial condition
    T =  initialize!(T, Tinitial)
    T = BC!(T, T_top = T_top, T_left = T_left, T_down = T_down, T_right = T_right)
    # main loop
    for t = 0:max_time
        Tnew = body!(T, Tnew, alpha = alpha, dt = dt, dx = dx, dy = dy)
        Tnew = BC!(Tnew, T_top = T_top, T_left = T_left, T_down = T_down, T_right = T_right)
        T = timeMarche!(T, Tnew)
    end
    return T
end
end # end of module FDM_D2Q4