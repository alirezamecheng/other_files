L = 1
N_f = 100
N_l = N_f
tau = 2.1
alpha = 5*10^-6

tstar = .1
dx_f = L/N_f
dt_f = (1/4)*(dx_f^2)/(alpha)

dt_l = 1
dx_l = 1

# tau = 2*alpha*dt_l / dx_l^2 + 0.5
alpha_l = ((dx_l)^2/(2*dt_l))*(tau - 0.5)
# LBMTimeSteps = tstar * N_l^2 * dx_l^2 / (alpha_l * dt_l)
LBMTimeSteps = 2 * tstar * N_l^2 /(tau - 0.5)
FDMTimeSteps = 4 * tstar * N_f^2

println("Time = ", tstar)
println("LBMTimeSteps = ",LBMTimeSteps)
println("FDMTimeSteps = ", FDMTimeSteps)
