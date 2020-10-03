export Transform

struct Transform{ℱ, ℬ}
    forward::ℱ
    backward::ℬ
end

function Transform(𝒢::FourierGrid)
    grid_size = length.(𝒢.grid)
    f = randn(grid_size...) .+ 0im
    FFTW.set_num_threads(Threads.nthreads())
    P = plan_fft(f) # typeof(iP) <: AbstractFFTs.ScaledPlan
    iP = plan_ifft(f)
    return Transform(P, iP)
end