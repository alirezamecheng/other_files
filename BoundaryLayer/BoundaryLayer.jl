#= 
        In the name of God.
    Boundary layer simulation in LBM.
    By: Alireza Ghavami nia. 
        
=# 



#/==============================================================================
#                                   Changes!
#/==============================================================================

#=
    * version 01 
        1- start of code.
        
        
=#


#/==============================================================================
#               Definition of constants and global variables
#/==============================================================================

# const Grid = 129
# const Nx = Grid - 1       # number of cells in the x-direction -1
# const Ny = Grid - 1       # number of cells in the y-direction -1
const Nx = 200       # number of cells in the x-direction -1
const Ny = 100       # number of cells in the y-direction -1
const Re = 100.0
MAX_ITTERATION_NUMBER = 10;
PrintPLTIntervals = 1


#/==============================================================================
#                             Choosing Method!
#/==============================================================================
 const MethodName = "BGK";                 # BGK Simple Boundary


println("Your Selected Method for solver is :",MethodName);
#/==============================================================================
#                       Printing And I/O Parameters
#/==============================================================================

const PrintDecimal = 1000000; # number of Decimals for file names.


#/==============================================================================
#                     Geometry VARIABLES and definition
#/==============================================================================
# WARNING: DO NOT Change ANY line below unless otherwise you are 200% sure what
# you are doing.

const Nx1 = Nx + 1;   # actual number of nodes in directions.
const Ny1 = Ny + 1;
const L = Ny + 1;     # width of the cavity
const Q = 9;          # number of discrete velocities
const rho0 = 1.0;     # initial density
const ux0 = 0.0;      # initial velocity component in x direction
const uy0 = 0.0;      # initial velocity component in y direction

#/==============================================================================
#                        BOUNDARY WALLS VELOCITIES
#/==============================================================================

const uw_top    = 0.0; # Lid U. Velocity or top X velocity
const vw_top    = 0.0; # Lid V. Velocity or top Y velocity

const uw_left   = 0.1; # Left wall U
const vw_left   = 0.0; # Left wall V

const uw_right  = 0.1; # Right wall U
const vw_right  = 0.0; # Right wall V

const uw_bottom = 0.0; # Bottom wall U
const vw_bottom = 0.0; # Bottom wall V

#/==============================================================================
#                             Methodes Constants
#/==============================================================================
#/=========================== public variables =================================

const cx = [1, 0, -1, 0, 1, -1, -1, 1, 0]; # lattice velocities changed 0-->9
const cy = [0, 1, 0, -1, 1, 1, -1, -1, 0]; # lattice velocities changed 0-->9

#/======================= BGK method constants =================================
const w  = [1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 4.0/9]; # lattice weights changed 0-->9






#/==============================================================================
#                             Used Packages
#/==============================================================================
# is needed just once 
# using Pkg
# Pkg.add("Dates")
# Pkg.add("Printf")



# alwayes is needed.
using Dates
using Printf









#/==============================================================================
#                           Results files Address
#/==============================================================================

const resultsAddress = "./Results/$(MethodName)_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))"; # for using in cluster.

mkpath(resultsAddress); # to creat specified path


println("Re = ",Re," Ux = ",uw_top," (Nx1 , Ny1) ==> (",Nx1,",",Ny1,")");
flush(stdout)











# Calculation the equilibrium distribution
function feq(RHO::Float64,U::Float64,V::Float64,k::Int64)  # RHO,U,V,k in this function are "Numbers".
    cu = cx[k]*U + cy[k]*V; # c K*u
    U2 = U*U + V*V;         # u*u
    return w[k]*RHO*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*U2);
end # end function [feq]


# Initialization
function Init_Eq!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64})
# چون یک بار اجرا می شود ارزش این که بهینه تعریف شود را ندارد
    @inbounds for c = 1:Nx1
      for r = 1:Ny1
        rho[r,c] = rho0;
        ux[r,c] = ux0;
        uy[r,c] = uy0;
       @simd for k = 1:Q
            f[r,c,k] = feq(rho[r,c], ux[r,c], uy[r,c], k); # recall feq function
        end #end for K
      end #end for r
    end #end for c
    return rho, ux, uy, f;
end #end function Init_Eq.



# BGK collision
function Coll_BGK!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64}, tau::Float64)
    #double FEQ
 @inbounds for k = 1:Q
        for c = 1:Nx1
           @simd for r = 1:Ny1
                FEQ = feq(rho[r,c], ux[r,c], uy[r,c], k); #EDF
                f_post[r,c,k] = f[r,c,k] - (f[r,c,k] - FEQ)/tau; #Post-collision DFs
            end # for r
        end # for c
    end #for k
    return f_post;
end #function [Coll_BGK]


# Streaming
function Streaming!(f::Array{Float64},f_post::Array{Float64})
  @inbounds for k = 1:Q
        for c = 1:Nx1
            for r = 1:Ny1
                cd = c - cx[k] # upwind node
                rd = r - cy[k]
                if cd>0 && rd>0 && cd<=Nx1 && rd<=Ny1
                    f[r,c,k] = f_post[rd,cd,k]; # streaming
                end # if
            end # for r
        end # for c
    end # for k
    return f;
end # function [Streaming]

# Bounce-back scheme
function Boundary_Conditions!(f::Array{Float64},f_post::Array{Float64},rho::Array{Float64})

     #  C = 1: left wall
     
#     @inbounds for r = 1:Ny1
#        f[r,1,1] = f_post[r,1,3];
#        f[r,1,5] = f_post[r,1,7]; # 1 singular node
#        f[r,1,8] = f_post[r,1,6]; # 1 singular node
#      end


    for r = 1:Ny1
        rhow = 1/(1-uw_left)*(f_post[r,1,9] + f_post[r,1,2] + f_post[r,1,4] + 2*(f_post[r,1,3] + f_post[r,1,6] + f_post[r,1,7]));
        f[r,1,1] = f_post[r,1,3] + 2/3*rhow*uw_left;
        f[r,1,5] = f_post[r,1,7] - 0.5*(f_post[r,1,2] - f_post[r,1,4]) + 1/6*rhow*uw_left + 0.5*rhow*vw_left;
        f[r,1,8] = f_post[r,1,6] + 0.5*(f_post[r,1,2] - f_post[r,1,4]) + 1/6*rhow*uw_left - 0.5*rhow*vw_left;
    end 
    
     # c = Nx: right wall
#     @inbounds for r = 1:Ny1
#        f[r,Nx1,3] = f_post[r,Nx1,1];
#        f[r,Nx1,7] = f_post[r,Nx1,5]; # 1 singular node
#        f[r,Nx1,6] = f_post[r,Nx1,8]; # 1 singular node
#      end

    for r = 1:Ny1
        rhow = 1/(1+uw_right)*(f_post[r,Nx1,9] + f_post[r,Nx1,2] + f_post[r,Nx1,4] + 2*(f_post[r,Nx1,1] + f_post[r,Nx1,5] + f_post[r,Nx1,8]));
        f[r,Nx1,3] = f_post[r,Nx1,1] - 2/3*rhow*uw_right;
        f[r,Nx1,7] = f_post[r,Nx1,5] + 0.5*(f_post[r,Nx1,2] - f_post[r,Nx1,4]) - 1/6*rhow*uw_right - 0.5*rhow*vw_right;
        f[r,Nx1,6] = f_post[r,Nx1,8] - 0.5*(f_post[r,Nx1,2] - f_post[r,Nx1,4]) - 1/6*rhow*uw_right + 0.5*rhow*vw_right;
    end 

     # r = Ny top plate
    @inbounds for c = 1:Nx1
         f[Ny1,c,4] = f_post[Ny1,c,2];
         f[Ny1,c,7] = f_post[Ny1,c,5]+6*rho[Ny1,c]*w[7]*cx[7]*uw_top;
         f[Ny1,c,8] = f_post[Ny1,c,6]+6*rho[Ny1,c]*w[8]*cx[8]*uw_top;
     end

     # r = 1: bottom plate
    @inbounds for c = 1:Nx1
       f[1,c,2] = f_post[1,c,4];
       f[1,c,5] = f_post[1,c,7];
       f[1,c,6] = f_post[1,c,8];
     end

     return f;
end # function [Boundary_Conditions!]


# Fluid variables (density and velocity)
function Den_vel!(rho::Array{Float64},ux::Array{Float64},uy::Array{Float64},f::Array{Float64})
  @inbounds for c = 1:Nx1
     @simd for r = 1:Ny1
            rho[r,c]=(f[r,c,9]+f[r,c,1]+f[r,c,2])+(f[r,c,3]+f[r,c,4]+f[r,c,5])+(f[r,c,6]+f[r,c,7]+f[r,c,8]);
            ux[r,c]=(f[r,c,1]+f[r,c,5]+f[r,c,8]-f[r,c,3]-f[r,c,6]-f[r,c,7])/rho[r,c];
            uy[r,c]=(f[r,c,5]+f[r,c,6]+f[r,c,2]-f[r,c,7]-f[r,c,8]-f[r,c,4])/rho[r,c];
        end # for r
    end # for c
    return rho , ux , uy;
end # function  [Den_vel!]

# function for Data output in plt format
function Data_output(LoopItterateCounter::Int64, rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64})
    FileNameAndAddress = @sprintf("%s/Results%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal)
    File = open(FileNameAndAddress,"w");
      println(File,"TITLE=\"",MethodName," Re = ",Re,"\"");
      # write(File,"SOLUTIONTIME = $(LoopItterateCounter)\n");
      # write(File,"T = \"Time=$(LoopItterateCounter)\"\n");
      write(File,"VARIABLES = X,Y,Rho,U,V\n");
      write(File,"ZONE\n");
      println(File,"I = ",Nx1,"\t","J = ",Ny1);
      write(File,"F=POINT\n");
    @inbounds  for r = 1:Ny1        # important note : This nested for loop should be like this
          for c = 1:Nx1
              # write(File,@sprintf("%d,%d,%lf,%lf,%lf\n",c,r,rho[r,c],ux[r,c],uy[r,c]));
              println(File,c,",",r,",",rho[r,c],",",ux[r,c],",",uy[r,c]);
          end
      end
    close(File);
    println("Results saved at : ",FileNameAndAddress);
    flush(stdout);
end



#///////////////////////////////////////////////////////////////////////////////
#
##                              Function main
#
#///////////////////////////////////////////////////////////////////////////////

function main()

    #/==============================================================================
    #                              Allocations
    #/==============================================================================

    f       = Array{Float64}(undef,Ny1,Nx1,Q); # array of the distribution functions (DFs)
    f_post  = Array{Float64}(undef,Ny1,Nx1,Q);

    rho     = Array{Float64}(undef,Ny1,Nx1);   # arrays of fluid density and velocity
    ux      = Array{Float64}(undef,Ny1,Nx1);
    uy      = Array{Float64}(undef,Ny1,Nx1);

    u0      = Array{Float64}(undef,Ny1,Nx1);   # for calculations in Err function
    v0      = Array{Float64}(undef,Ny1,Nx1);


    #/================= SRT Special Variables and allocations ======================

    tau::Float64     = 3*L*uw_top/Re + 0.5;      # relaxation time for BGK



    #/==============================================================================
    #                   Initialization of rho, ux, uy and f
    #/==============================================================================

    rho, ux, uy, f = Init_Eq!(rho, ux, uy, f);

    

    for LoopItterateCounter=1:MAX_ITTERATION_NUMBER
            f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);               # BGK collision
            f = Streaming!(f,f_post);                                      # Streaming for BGK
            f = Boundary_Conditions!(f,f_post,rho);                        # Boundary Conditions
            rho, ux, uy = Den_vel!(rho,ux,uy,f);                           # Fluid variables calculation for all models
        
    
    if LoopItterateCounter%PrintPLTIntervals ==0
        Data_output(LoopItterateCounter,rho,ux,uy);       # Output last simulation data
        println("LoopItterateCounter: ",LoopItterateCounter)

    end # end of if [PrintPLTIntervals]
    
    
    
    end # end main loop 

    
end # end main function
    

@time main()
