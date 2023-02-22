using Plots

include("fgn.jl")

function fOU(m, α, v, H, N, T; X₀=0, npath=1, seed=nothing)
    """
    fraction ornstein-uhlenbeck base log volatility
    """
    # TODO, make a params struct or make params a named tuple
        # if struct move bounds checking there
    @assert α>0 "α must be >0"
    @assert v>0 "v must be >0"

    Wᴴ = fractional_gaussian_noise(H,N,T; npath=npath,seed=seed)
    δt = T/N
    Xₜ = zeros(Float64, (N-1,npath))
    if npath>1
        Xₜ[1,:] .= X₀
    else
        Xₜ[1] = X₀
    end
    for j in 1:npath
        for i in 2:N-1
            Xₜ[i,j] = α*(m - Xₜ[i-1,j])*δt + v*(Wᴴ[i,j] - Wᴴ[i-1,j])
        end
    end
    return Xₜ
end

function fVol()
    
end