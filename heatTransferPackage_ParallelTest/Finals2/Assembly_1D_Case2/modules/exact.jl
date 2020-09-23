@info "Loading Plots"
using Plots
@info "Loaded."
function cn(n,L)
    return (2/L)*( (n * pi) / L * (-1)^n + (n^2 * pi ^2)/(L^2) * (-1)^n) - (n^2* pi^2)/(L^2)
end

function main()

    L = 0.05
    h = 0.00125
    max_time = 2
    landa = 35
    c = 4.875 * 10^6
    Tw = 0
    Te = 100
    Tinitial = 0

    dx = h
    alpha = landa / c
    Nx = Int64(L / dx + 1)
    x = 0:dx:L
    Ts1 = Tw
    Ts2 = Te
    t = max_time
    T = Array{Float64,1}(undef,size(x,1))
    T .= 0
    part = 0
    @info "Solving"
    for i = 1:size(x,1)
        for n = 1:1
            mylanda = n*pi/L
            part = part + cn(n,L) * sin(pi*n/L) * exp(-alpha * mylanda^2 * t) 
        end
        T[i] = part + (Ts2 - Ts1)/L * x[i] + Ts1
        "T = $(T)"
    end
    @info "Ploting"
    plot(x,T)
end

main()