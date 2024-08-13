include("fgn.jl")

struct σParams(s, x, a) 
    ### toying around with this idea
    # can use a struct to put bounds on params probably
    σ₀ = s
    ξ = x
    α = a
end

function sigma_t(fgn::T, t, σ₀, ξ; α=1) where T <: AbstractArray
    r(s) = s^(2H)
    N = size(fgn)[1]
    δt = t/N
    σₜ = [σ₀, zeros(N)...]
	for i in 1:N
		rₜ = r(δt * i)
		Bᴴ = fgn[i]
		σₜ[i+1] = σ₀ * exp(ξ*Bᴴ - 0.5α * ξ^2 * rₜ)
	end
	return σₜ
end

function sigma_t(H, N, T, σ₀, ξ; α=1, nsim=1, seed=1234)
    fgn = fractional_gaussian_noise(H, N, t; nsim=nsim, seed=seed)
    δt = T/N
    σₜ = [σ₀, zeros(N)...]
	for i in 1:N
		rₜ = r(δt * i)
		Bᴴ = fgn[i]
		σₜ[i+1] = σ₀ * exp(ξ*Bᴴ - 0.5α * ξ^2 * rₜ)
	end
	return σₜ
end

begin
    H = 0.1
	q = 10
	n = 2^q + 1
	t = 1
    fgn = fractional_gaussian_noise(H,n,t)

	σ₀ = 0.15
	ξ = 0.8
	α = 1
    σ = sigma_t(fgn, 1, σ₀, ξ; α=α)
end;