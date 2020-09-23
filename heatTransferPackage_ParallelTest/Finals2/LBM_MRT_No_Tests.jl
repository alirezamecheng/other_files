
# Constants
const N = 10

const Q = 9
const W = [1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36, 4/9]
const Nx = N
const Ny = N

const cx = [1, 0, -1, 0, 1, -1, -1, 1, 0] # lattice velocities changed 0-->9
const cy = [0, 1, 0, -1, 1, 1, -1, -1, 0] # lattice velocities changed 0-->9

function Collision!(g::Array{Float64}, gpost::Array{Float64}, S::Array{Float64}, T::Array{Float64})
    S0 = S[9]
    S1 = S[1]
    S2 = S[2]
    S3 = S[3]
    S4 = S[4]
    S5 = S[5]
    S6 = S[6]
    S7 = S[7]
    S8 = S[8]

    for j = 1:Ny
        for i = 1:Nx

            Temp = T[i,j];

            g0 = g[i,j,9];
            g1 = g[i,j,1];
            g2 = g[i,j,2];
            g3 = g[i,j,3];
            g4 = g[i,j,4];
            g5 = g[i,j,5];
            g6 = g[i,j,6];
            g7 = g[i,j,7];
            g8 = g[i,j,8];
            M0=(     g0)+(     g1)+(     g2)+(     g3)+(     g4)+(     g5)+(     g6)+(     g7)+(     g8);
            M1=(-4.0*g0)+(-1.0*g1)+(-1.0*g2)+(-1.0*g3)+(-1.0*g4)+( 2.0*g5)+( 2.0*g6)+( 2.0*g7)+( 2.0*g8);
            M2=( 4.0*g0)+(-2.0*g1)+(-2.0*g2)+(-2.0*g3)+(-2.0*g4)+(     g5)+(     g6)+(     g7)+(     g8);
            M3=( 0.0*g0)+(     g1)+( 0.0*g2)+(-1.0*g3)+( 0.0*g4)+(     g5)+(-1.0*g6)+(-1.0*g7)+(     g8);
            M4=( 0.0*g0)+(-2.0*g1)+( 0.0*g2)+( 2.0*g3)+( 0.0*g4)+(     g5)+(-1.0*g6)+(-1.0*g7)+(     g8);
            M5=( 0.0*g0)+( 0.0*g1)+(     g2)+( 0.0*g3)+(-1.0*g4)+(     g5)+(     g6)+(-1.0*g7)+(-1.0*g8);
            M6=( 0.0*g0)+( 0.0*g1)+(-2.0*g2)+( 0.0*g3)+( 2.0*g4)+(     g5)+(     g6)+(-1.0*g7)+(-1.0*g8);
            M7=( 0.0*g0)+(     g1)+(-1.0*g2)+(     g3)+(-1.0*g4)+( 0.0*g5)+( 0.0*g6)+( 0.0*g7)+( 0.0*g8);
            M8=( 0.0*g0)+( 0.0*g1)+( 0.0*g2)+( 0.0*g3)+( 0.0*g4)+(     g5)+(-1.0*g6)+(     g7)+(-1.0*g8);

            Meq0 = geq(Temp,9)
            Meq1 = geq(Temp,1)
            Meq2 = geq(Temp,2)
            Meq3 = geq(Temp,3)
            Meq4 = geq(Temp,4)
            Meq5 = geq(Temp,5)
            Meq6 = geq(Temp,6)
            Meq7 = geq(Temp,7)
            Meq8 = geq(Temp,8)

            # M0_d=M0+(S0*(Meq0-M0));
            # M1_d=M1+(S1*(Meq1-M1));
            # M2_d=M2+(S2*(Meq2-M2));
            # M3_d=M3+(S3*(Meq3-M3));
            # M4_d=M4+(S4*(Meq4-M4));
            # M5_d=M5+(S5*(Meq5-M5));
            # M6_d=M6+(S6*(Meq6-M6));
            # M7_d=M7+(S7*(Meq7-M7));
            # M8_d=M8+(S8*(Meq8-M8));

            M0_d=(S0*(Meq0-M0));
            M1_d=(S1*(Meq1-M1));
            M2_d=(S2*(Meq2-M2));
            M3_d=(S3*(Meq3-M3));
            M4_d=(S4*(Meq4-M4));
            M5_d=(S5*(Meq5-M5));
            M6_d=(S6*(Meq6-M6));
            M7_d=(S7*(Meq7-M7));
            M8_d=(S8*(Meq8-M8));

            gpost[i,j,9] =g0 + ((1.0/9.0)*M0_d)+((-4.0/36.0)*M1_d)+(( 4.0/36.0)*M2_d)+(( 0.0/6.0)*M3_d)+(( 0.0/12.0)*M4_d)+(( 0.0/6.0)*M5_d)+(( 0.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);
            gpost[i,j,1] =g1 + ((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+(( 1.0/6.0)*M3_d)+((-2.0/12.0)*M4_d)+(( 0.0/6.0)*M5_d)+(( 0.0/12.0)*M6_d)+(( 1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);
            gpost[i,j,2] =g2 + ((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+(( 0.0/6.0)*M3_d)+(( 0.0/12.0)*M4_d)+(( 1.0/6.0)*M5_d)+((-2.0/12.0)*M6_d)+((-1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);
            gpost[i,j,3] =g3 + ((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+((-1.0/6.0)*M3_d)+(( 2.0/12.0)*M4_d)+(( 0.0/6.0)*M5_d)+(( 0.0/12.0)*M6_d)+(( 1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);
            gpost[i,j,4] =g4 + ((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+(( 0.0/6.0)*M3_d)+(( 0.0/12.0)*M4_d)+((-1.0/6.0)*M5_d)+(( 2.0/12.0)*M6_d)+((-1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);
            gpost[i,j,5] =g5 + ((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+(( 1.0/6.0)*M3_d)+(( 1.0/12.0)*M4_d)+(( 1.0/6.0)*M5_d)+(( 1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+(( 1.0/4.0)*M8_d);
            gpost[i,j,6] =g6 + ((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+((-1.0/6.0)*M3_d)+((-1.0/12.0)*M4_d)+(( 1.0/6.0)*M5_d)+(( 1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+((-1.0/4.0)*M8_d);
            gpost[i,j,7] =g7 + ((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+((-1.0/6.0)*M3_d)+((-1.0/12.0)*M4_d)+((-1.0/6.0)*M5_d)+((-1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+(( 1.0/4.0)*M8_d);
            gpost[i,j,8] =g8 + ((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+(( 1.0/6.0)*M3_d)+(( 1.0/12.0)*M4_d)+((-1.0/6.0)*M5_d)+((-1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+((-1.0/4.0)*M8_d);
    end #end [for i]
end #end [for j]

return gpost
end # function [Coll_MRT2]

function geq(T::Float64, q::Int)
    return W[q] * T
end

# function Streaming!(g::Array{Float64},gpost::Array{Float64})
#     for k = 1:Q
#         for j = 1:Ny
#             for i = 1:Nx
#                 id = i - cx[k] # upwind node
#                 jd = j - cy[k]
#                 if id>0 && jd>0 && id<=Nx && jd<=Ny
#                     g[i,j,k] = gpost[id,jd,k]; # streaming
#                 end # if
#             end # for i
#         end # for j
#     end # for k
#     return g;
# end # function [Streaming]
function Streaming!(g::Array{Float64},gpost::Array{Float64})
    #=
    lattice configuration 
       6    2    5
        \   ^   /
         \  |  /
          \ | /
    3 <-----0-----> 1
           /|\
          / | \
         /  |  \
       7    4   8
    =#
    k = 9 
    for j = 1:Ny
        for i = 1:Nx
            g[i,j,k] = gpost[i,j,k]
        end
    end
    k = 1 
    for j = 1:Ny
        for i = Nx-1:-1:1
            g[i+1,j,k] = gpost[i,j,k]
        end
    end
    k = 2 
    for j = Ny-1:-1:1
        for i = 1:1:Nx
            g[i,j+1,k] = gpost[i,j,k]
        end
    end
    k = 3 
    for j = 1:Ny
        for i = 2:Nx
            g[i-1,j,k] = gpost[i,j,k]
        end
    end
    k = 4 
    for j = 2:Ny
        for i = 1:Nx
            g[i,j-1,k] = gpost[i,j,k]
        end
    end
    k = 5 
    for j = Ny-1:-1:1
        for i = Nx-1:-1:1
            g[i+1,j+1,k] = gpost[i,j,k]
        end
    end
    k = 6 
    for j = Ny-1:-1:1
        for i = 2:Nx
            g[i-1,j+1,k] = gpost[i,j,k]
        end
    end
    k = 7 
    for j = 2:Ny
        for i = 2:Nx
            g[i-1,j-1,k] = gpost[i,j,k]
        end
    end
    k = 8 
    for j = 2:Ny
        for i = Nx-1:-1:1
            g[i+1,j-1,k] = gpost[i,j,k]
        end
    end
    return g
end

function BC!(g::Array{Float64}, gpost::Array{Float64}; T_top::Real, T_right::Real, T_down::Real, T_left::Real)
    # top BC
    j = Ny
    for i = 1:Nx
        g[i,j,1] = -gpost[i,j,3] + (2/9)*T_top
        g[i,j,2] = -gpost[i,j,4] + (2/9)*T_top
        g[i,j,3] = -gpost[i,j,1] + (2/9)*T_top
        g[i,j,4] = -gpost[i,j,2] + (2/9)*T_top
        g[i,j,5] = -gpost[i,j,7] + (1/18)*T_top
        g[i,j,6] = -gpost[i,j,8] + (1/18)*T_top
        g[i,j,7] = -gpost[i,j,5] + (1/18)*T_top
        g[i,j,8] = -gpost[i,j,6] + (1/18)*T_top
    end

    # right BC
    i = Nx
    for j = 1:Ny
        g[i,j,1] = -gpost[i,j,3] + (2/9)*T_right
        g[i,j,2] = -gpost[i,j,4] + (2/9)*T_right
        g[i,j,3] = -gpost[i,j,1] + (2/9)*T_right
        g[i,j,4] = -gpost[i,j,2] + (2/9)*T_right
        g[i,j,5] = -gpost[i,j,7] + (1/18)*T_right
        g[i,j,6] = -gpost[i,j,8] + (1/18)*T_right
        g[i,j,7] = -gpost[i,j,5] + (1/18)*T_right
        g[i,j,8] = -gpost[i,j,6] + (1/18)*T_right
    end

    # down BC
    j = 1
    for i = 1:Nx
        g[i,j,1] = -gpost[i,j,3] + (2/9)*T_down
        g[i,j,2] = -gpost[i,j,4] + (2/9)*T_down
        g[i,j,3] = -gpost[i,j,1] + (2/9)*T_down
        g[i,j,4] = -gpost[i,j,2] + (2/9)*T_down
        g[i,j,5] = -gpost[i,j,7] + (1/18)*T_down
        g[i,j,6] = -gpost[i,j,8] + (1/18)*T_down
        g[i,j,7] = -gpost[i,j,5] + (1/18)*T_down
        g[i,j,8] = -gpost[i,j,6] + (1/18)*T_down
    end

    # left BC
    i = 1
    for j = 1:Ny
        g[i,j,1] = -gpost[i,j,3] + (2/9)*T_left
        g[i,j,2] = -gpost[i,j,4] + (2/9)*T_left
        g[i,j,3] = -gpost[i,j,1] + (2/9)*T_left
        g[i,j,4] = -gpost[i,j,2] + (2/9)*T_left
        g[i,j,5] = -gpost[i,j,7] + (1/18)*T_left
        g[i,j,6] = -gpost[i,j,8] + (1/18)*T_left
        g[i,j,7] = -gpost[i,j,5] + (1/18)*T_left
        g[i,j,8] = -gpost[i,j,6] + (1/18)*T_left
    end
    return g
end

function Tcalc!(T::Array{Float64}, g::Array{Float64})
    for j = 1:Ny
        for i = 1:Nx 
            T[i,j] = g[i,j,1] + g[i,j,2] + g[i,j,3] + g[i,j,4] + g[i,j,5] + g[i,j,6] + g[i,j,7] + g[i,j,8] + g[i,j,9]
        end
    end
    return T
end

function initialCondition!(T::Array{Float64}, g::Array{Float64})
    for q = 1:Q
        for j = 1:Ny
            for i = 1:Nx
                g[i,j,q] = geq(T[i,j], q)
            end
        end
    end
    return g
end

function main()
    # problem Spec.
    # Max_Time_steps = 100
    Max_Time = 2000 #[s]
    dt = 1.0
    T_t = 10
    T_r = 10
    T_d = 10
    T_l = 10
    D = 0.01 # diffusion Coefficient 

    L = 0.1 # [m]
    # dx = L/(Nx-1) 
    dx = 1.0

    # Allocations
    T = Array{Float64, 2}(undef, N, N)
    g = Array{Float64, 3}(undef, N, N, Q)
    gpost = Array{Float64, 3}(undef, N, N, Q)

    # Coefficients Calculation 
    S = Array{Float64}(undef,Q)
    m = Array{Float64}(undef,Q)

    # S[1] = 1
    # S[2] = 1
    # S[3] = 1 # convective Coefficient
    # S[4] = 1 # convective Coefficient
    # S[5] = 1 # convective Coefficient
    # S[6] = 1 # convective Coefficient
    # S[7] = 1
    # S[8] = 1
    # S[9] = 1 # the relaxation coefficient for the conserved scalar variable
    S .= 0.6
    # @show S .= 0.5 + dt*D*3/(dx^2)

    # initial condition 
    @show T .= 0
    g = initialCondition!(T, g)
    gpost .= g
    g = BC!(g, gpost, T_top = T_t, T_right = T_r, T_down = T_d, T_left = T_l)

    # main Loop 
    for t = 0:dt:Max_Time
        gpost = Collision!(g, gpost, S, T)
        g = Streaming!(g, gpost)
        g = BC!(g, gpost, T_top = T_t, T_right = T_r, T_down = T_d, T_left = T_l)
        T = Tcalc!(T, g)
    end # end main Loop [t]

    T
end

a = main()
