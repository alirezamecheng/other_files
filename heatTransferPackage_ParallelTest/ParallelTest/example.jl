@info("Compiling packages")
using BenchmarkTools, Statistics, DataFrames, CSV, Distributed, LinearAlgebra, SharedArrays
@info("Number of Threads = ",Threads.nthreads())

function naive_gemv(alpha, A, x, beta, y)
    (alpha * A * x) + (beta * y)
end

function sequential_gemv!(alpha::Float64,
                          A::Array{Float64, 2},
                          x::Array{Float64, 1},
                          beta::Float64,
                          y::Array{Float64, 1},
                          output::Array{Float64, 1})

    @inbounds for i in 1:size(A, 1)
        @inbounds for j in 1:size(x, 1)
            output[i] += A[i, j] * x[j]
        end
        output[i] *= alpha
        output[i] += beta * y[i]
    end

    output
end

function threads_gemv!(alpha::Float64,
                       A::Array{Float64, 2},
                       x::Array{Float64, 1},
                       beta::Float64,
                       y::Array{Float64, 1},
                       output::Array{Float64, 1})

    @inbounds Threads.@threads for i in 1:size(A, 1)
        @inbounds Threads.@threads for j in 1:size(x, 1)
            output[i] += A[i, j] * x[j]
        end
        output[i] *= alpha
        output[i] += beta * y[i]
    end

    output
end

function processes_gemv!(alpha::Float64,
                         A::Array{Float64, 2},
                         x::Array{Float64, 1},
                         beta::Float64,
                         y::Array{Float64, 1},
                         output::SharedArray{Float64, 1})

    @inbounds @sync @distributed for i in 1:size(A, 1)
        @inbounds for j in 1:size(x, 1)
            output[i] += A[i, j] * x[j]
        end
        output[i] *= alpha
        output[i] += beta * y[i]
    end

    output
end

function benchmark_gemv()
    samples = 8

    alpha = 0.5
    beta = 0.2

    starting_m = 2 ^ 8
    max_m = 2 ^ 15

    starting_n = 2 ^ 9
    max_n = 2 ^ 9

    m = starting_m
    n = starting_n

    functions_number = 5
    allowed_duration = 30

    results = Array{Float64, 2}(undef,
                                functions_number * (Int(log2(max_m) - log2(starting_m) + 1)),
                                samples + 5)
    fill!(results, 0.0)

    y = rand(n)
    x = rand(m)
    A = rand(n, m)
    output = similar(y)
    fill!(output, 0.0)

    shared_output = SharedArray(output)

    BLAS.set_num_threads(Threads.nthreads())

    if nprocs() <= 1
        addprocs(5)
    end

    println("Setting BLAS Threads to $(Threads.nthreads())")
    println("Running Threads: $(Threads.nthreads())")
    println("Running processes: $(nprocs())")

    println("Warmup runs:")
    println("Running naive_gemv")
    @benchmark naive_gemv($alpha, $A, $x, $beta, $y) samples=1 gcsample=true seconds=5

    println("Running sequential_gemv!")
    @benchmark sequential_gemv!($alpha, $A, $x, $beta, $y, $output) samples=1 gcsample=true seconds=5

    println("Running threads_gemv!")
    @benchmark threads_gemv!($alpha, $A, $x, $beta, $y, $output) samples=1 gcsample=true seconds=5

    # println("Running processes_gemv!")
    # @benchmark processes_gemv!($alpha, $A, $x, $beta, $y, $shared_output) samples=1 gcsample=true seconds=5

    println("Running BLAS.gemv!")
    @benchmark BLAS.gemv!('N', $alpha, $A, $x, $beta, $y) samples=samples gcsample=true seconds=5

    println("Starting benchmark")
    for i = 1:functions_number:(functions_number * (Int(log2(max_m) - log2(starting_m) + 1)))
        println("Iteration $i")

        y = rand(n)
        x = rand(m)
        A = rand(n, m)
        output = similar(y)
        fill!(output, 0.0)

        allowed_seconds = (1 + div(i, functions_number)) * allowed_duration

        println("Running naive_gemv")
        result = @benchmark naive_gemv($alpha, $A, $x, $beta, $y) samples=samples gcsample=true seconds=allowed_seconds
        results[i, :] = vcat(result.times, alpha, beta, m, n, 1)

        println("Running sequential_gemv!")
        result = @benchmark sequential_gemv!($alpha, $A, $x, $beta, $y, $output) samples=samples gcsample=true seconds=allowed_seconds
        results[i + 1, :] = vcat(result.times, alpha, beta, m, n, 2)

        fill!(output, 0.0)

        println("Running threads_gemv!")
        result = @benchmark threads_gemv!($alpha, $A, $x, $beta, $y, $output) samples=samples gcsample=true seconds=allowed_seconds
        results[i + 2, :] = vcat(result.times, alpha, beta, m, n, 3)

        fill!(output, 0.0)
        shared_output = SharedArray(output)

        println("Running processes_gemv!")
        result = @benchmark processes_gemv!($alpha, $A, $x, $beta, $y, $shared_output) samples=samples gcsample=true seconds=allowed_seconds
        results[i + 3, :] = vcat(result.times, alpha, beta, m, n, 4)

        println("Running BLAS.gemv!")
        result = @benchmark BLAS.gemv!('N', $alpha, $A, $x, $beta, $y) samples=samples gcsample=true seconds=allowed_seconds
        results[i + 4, :] = vcat(result.times, alpha, beta, m, n, 5)

        m *= 2
        # We're keeping n fixed
        #n *= 2
    end

    results
end
benchmark_gemv()