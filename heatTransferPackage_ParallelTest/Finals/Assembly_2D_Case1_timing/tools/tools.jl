function saveTP(T,name)
    Nx = size(T,1)
    Ny = size(T,2)
    mkpath("results")
    io = open("results/$(name).plt","w")
    println(io,"TITLE = \"$(name)\"")
    println(io,"VARIABLES = X,Y,T")
    println(io,"ZONE")
    println(io,"I = ",Nx,"\t","J = ",Ny)
    println(io,"F=POINT")
    for j = 1:Ny
        for i = 1:Nx
            println(io,i,",",j,",",T[i,j])
        end 
    end
    close(io)
    println("Results are saved at: Result_$(name)")
end
function saveCenterLines(T,name; Orientation = :x)
    mkpath("results")
    if Orientation == :x
        N = size(T,1)
        orient = "X"
        TT = T[:, Int64(floor(size(T,2)/2))]
    else Orientation == :y
        N = size(T,2)
        orient = "Y"
        TT = T[Int64(floor(size(T,1)/2)),:]
    end
    io = open("results/cneter_$(orient)_$(name).plt","w")
    println(io,"TITLE = \"$(name)\"")
    println(io,"VARIABLES = $(orient),T")
    println(io,"ZONE")
    println(io,"I = ",N)
    println(io,"F=POINT")
    for i = 1:N
            println(io,i,",",TT[i])
    end 
    close(io)
    println("Results are saved at: Result_$(orient)_$(name)")
end
