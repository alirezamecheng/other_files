# D2Q4 from [A.A.Mohammad]
# const w = [0.5 0,5] are hard coded in the lattice
# There is no Allocation for `fpost` in this code

module LBM_BGK_D2Q4
export solve, LBM_BGK_D2Q4_Problem_maker
struct LBM_BGK_D2Q4_Problem
    N::Int64
    tau::Float64
    max_time::Int64
    Tinitial::Float64
    T_top::Float64
    T_left::Float64
    T_down::Float64
    T_right::Float64
    dx::Float64
    dt::Float64
end

function LBM_BGK_D2Q4_Problem_maker(;N::Int64,
                                tau::Number,
                                max_time::Int64,
                                Tinitial::Number,
                                T_top::Number,
                                T_left::Number,
                                T_down::Number,
                                T_right::Number,
                                dx = 1.0,
                                dt = 1.0)
    P = LBM_BGK_D2Q4_Problem(Int64(N),
                        Float64(tau),
                        Int64(max_time),
                        Float64(Tinitial),
                        Float64(T_top),
                        Float64(T_left),
                        Float64(T_down),
                        Float64(T_right),
                        Float64(dx),
                        Float64(dt))

    return P
end

function feq(T::Float64)
    return 0.25*T 
end
function collision!(f::Array{Float64}, T::Array{Float64}; dt::Float64, tau::Float64)
    for q = 1:4
        for j = 1:size(f,2)
            for i = 1:size(f,1)
                f[i,j,q] = f[i,j,q]*(1-(dt/tau)) + (dt/tau)*feq(T[i,j])
            end
        end
    end
    return f
end
function streaming!(f::Array{Float64})
    Nx = size(f,1)
    Ny = size(f,2)
    q = 1
    for j = 1:Ny 
        for i = (Nx-1):-1:1
            f[i+1,j,q] = f[i,j,q]
        end
    end
    q = 2
    for j = 1:Ny
        for i = 2:Nx
            f[i-1,j,q] = f[i,j,q]
        end
    end
    q = 3
    for j = (Ny-1):-1:1
        for i = 1:Nx
            f[i,j+1,q] = f[i,j,q]
        end
    end
    q = 4
    for j = 2:Ny
        for i = 1:Nx
            f[i,j-1,q] = f[i,j,q]
        end
    end
    return f    
end
function BC!(f::Array{Float64}; T_top::Float64, T_left::Float64, T_down::Float64, T_right::Float64)
    Nx = size(f,1)
    Ny = size(f,2)
    #=
    Lattice Orientation
            3
            ^
            |
            |
    2 <-----+-----> 1
            |
            |
            4
    =#
    # left  wall
    i = 1
    for j = 1:Ny
        f[i,j,1] = 0.5*T_left - f[i,j,2]
    end
    # right wall
    i = Nx
    for j = 1:Ny
        f[i,j,2] = 0.5*T_right - f[i,j,1]
    end
    # Top wall 
    j = Ny
    for i = 1:Nx
        f[i,j,4] = 0.5*T_top - f[i,j,3]
    end
    # bottom wall
    j = 1
    for i  = 1:Nx
        f[i,j,3] = 0.5*T_down - f[i,j,4]
    end
    return f
end
function initialize!(f::Array{Float64},Tinitial::Float64)
    for q = 1:4
        for j = 1:size(f,2)
            for i = 1:size(f,1)
                f[i,j,q] = feq(Tinitial)
            end
        end
    end
    return f
end
function Tcalc!(T::Array{Float64}, f::Array{Float64})
    for j = 1:size(f,2)
        for i = 1:size(f,1)
            T[i,j] = f[i,j,1] + f[i,j,2] + f[i,j,3] + f[i,j,4]
        end
    end
    return T
end
function solve(P::LBM_BGK_D2Q4_Problem)

    # parameters
    N = P.N
    tau = P.tau
    Nx = P.N
    Ny = P.N

    Tinitial = P.Tinitial
    T_top = P.T_top
    T_left = P.T_left
    T_down = P.T_down
    T_right = P.T_right

    # calculation of other parameters
    dt = P.dx 
    dx = P.dt
    # tau = 2 * alpha * dt / (dx^2) + 0.5
    alpha = dx^2/(2*dt^2) *(tau - 0.5)

    maxIteration = P.max_time

    # Allocations
    T = Array{Float64, 2}(undef, Nx, Ny)
    f = Array{Float64, 3}(undef, Nx, Ny, 4)

    # initialize
    T .= Tinitial
    f = initialize!(f,Float64(Tinitial))

    for t = 1:maxIteration
        f = collision!(f, T, dt = dt, tau = tau)
        f = streaming!(f)
        f = BC!(f, T_top = T_top, T_left = T_left, T_down = T_down, T_right = T_right)
        T = Tcalc!(T, f)
    end
    return T
end
end # end of module [BGK_LBM_D2Q4]
