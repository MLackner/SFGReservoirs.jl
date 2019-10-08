# TODO: not sure about phase contribution
"""
"""
function sfspectrum(x, A, ω, Γ; nr=0.0, offset=0.0, nrphase=0.0)
    y = fill(0.0 + 0.0im, length(x))
    yr = Array{Float64}(undef, length(x))
    for i = 1:size(A,1)
        @. y += A[i] / (x - ω[i] - 1im * Γ[i]) #* exp(1im * φ[i])
    end
    @. yr = abs2(nr * exp(1im * nrphase) + y) + offset
end
