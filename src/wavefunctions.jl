using LinearAlgebra, Random

include("../../Constants/src/Short.jl")

@doc raw"""
the free localized wave function is given as (1/√(2πħ))*∫ exp(i(p.x - Et)/ħ) dp
this example only takes into account one space dimension

# input
x::Float64

 ::Array{Float64,1} - the spacial values
t::Float64

 ::Array{Float64,1} - the time values

p::Float64

 ::Array{Float64,1} - the momentum values

m::Float64 - the mass of the particle

free::String - are we dealing with a free particle or not

localized::String - are we dealing with a localized particle or not

len::Int - the number of values we give to the integral over p

# output
ψ::Array{Complex,3} - x is horizontal

                      p is vertical

                      t is depth

ψ::Array{Complex,3} - x is vertical

                    t is horizontal

# moral
we remove the dependence on p if the particle is localized, no longer requiring it as a variable

now there is also no mechanism to get ϕ(p)
"""
function ψ(x::Float64,t::Float64=1.; p::Float64=x, m::Float64=1., free::String="y", localized::String="y",
            len::Int=1000, L::Float64=0.)
    rng = Random.seed!(42)
    if free == "y"
        E = sqrt(dot(p.*c,p.*c) + (m*c^2)^2)
        if localized == "n"
            ψ = [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:length(p), j in 1:length(x), t in t]

        elseif localized == "y"
            if L == 0.
                p = range(-abs(100*x), 100*x, length=len)
            else
                p = range(0, step=2π/L, length=len)
            end
            ϕ = rand(rng,len)
            N = 1/(2π*ħ)
            ψ = [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:len, j in 1:length(x), t in t]
            return reshape((1/sqrt(2π*ħ))*(1/len)*sum(ϕ.*ψ, dims=1),:,length(t))
        end
    end
end
function ψ(x::Array{Float64,1},t::Array{Float64,1}=[1.]; p::Array{Float64,1}=x, m::Float64=1., free::String="y",
            localized::String="y", len::Int=1000, L::Float64=0.)
    rng = Random.seed!(42)
    if free == "y"
        E = sqrt(dot(p.*c,p.*c) + (m*c^2)^2)
        if localized == "n"
            ψ = [exp.(im*(dot(p[i],x[j]) .- E*t) / ħ) for i in 1:length(p), j in 1:length(x), t in t]
            return ψ
        elseif localized == "y"
            if L == 0.
                p = collect(range(-abs.(100*maximum(x)), 100*maximum(x), length=len))
            else
                p = collect(range(0, step=2π/L, length=len))
            end
            ϕ = rand(rng,len)
            N = 1/(2π*ħ)
            ψ = collect(exp.(im*(dot.(p[i],x[j]) .- E*t) / ħ) for i in 1:len, j in 1:length(x), t in t[:])
            return reshape((1/sqrt(2π*ħ))*(1/len)*sum(ϕ.*ψ, dims=1),:,length(t))
        end
    end
end

@doc raw"""
the probability amplitude is given as the modulus squared of the wavefunction

# input
x::Float64
 ::Array{Float64,1} - the spacial values
t::Float64
 ::Array{Float64,1} - the time values
p::Float64
 ::Array{Float64,1} - the momentum values
m::Float64 - the mass of the particle
free::String - are we dealing with a free particle or not
localized::String - are we dealing with a localized particle or not

# output
P::Array{Complex,3} - x is horizontal
                      p is vertical
                      t is depth

# moral
the probability is 1 everywhere with the old wavefunction defn

# TODO
     - have to consider in detail if element wise multiplication is really what i want here, we have a dot product like operation
     instead of a cummulative modulus operation of probabilty
     - moreover, we could instead dictate the bilinear, make the denity matrix and take the sum of the diagonal
"""
function P(x::Float64,t::Float64=1.; p::Float64=x, m::Float64=1., free::String="y", localized::String="y",
        len::Int=1000, fn::String="ψ")
    if fn == "ψ"
        return real(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len).*
            conj(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len)))
    elseif fn == "ϕ"
        return real(ϕ(x; m=m, len=len).*conj(ϕ(x; m=m, len=len)))
    end
end
function P(x::Array{Float64,1},t::Array{Float64,1}=[1.]; p::Array{Float64,1}=x, m::Float64=1., free::String="y",
            localized::String="y", len::Int=1000, fn::String="ψ")
    if fn == "ψ"
        return real(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len).*
            conj(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len)))
    elseif fn == "ϕ"
        return real(ϕ(x; m=m, len=len).*conj(ϕ(x; m=m, len=len)))
    end
end

@doc raw"""
    the momentum space wavefunction, or more exactly the Fourier transform of the position wavefunction

    # input
    p::Float64
    ::Array{Float64,1} - the momentum values
    m::Float64 - the mass of the particle
    len::Int - the number of values we give to the integral over x

    # output
    ϕ::Array{Complex,2} - x is horizontal
"""
function ϕ(p::Float64; m::Float64=1., len::Int=1000)
    x = range(-abs(100*p), 100*p, length=len)
    ψ = rand(len)[1,:]
    return (1/sqrt(2π*ħ))*(1/len)*sum(ψ .* [exp(-im*(dot(p[i],x[j])) / ħ)
            for i in 1:length(p), j in 1:len], dims=2)
end
function ϕ(p::Array{Float64,1}; m::Float64=1., len::Int=1000)
    x = range(-abs(100*minimum(p)), 100*maximum(p), length=len)
    ψ = rand(len)[1,:]
    return (1/sqrt(2π*ħ))*(1/len)*sum(ψ .* [exp(-im*(dot(p[i],x[j])) / ħ)
            for i in 1:length(p), j in 1:len], dims=2)
end

@doc raw"""

"""
function bilinears(x::Float64,t::Float64=1.; p::Float64=x, m::Float64=1., free::String="y", localized::String="y",
        len::Int=1000, fn::String="ψ")
    if fn == "ψ"
        return ψ(x,t; p=p, m=m, free=free, localized=localized, len=len).*
            conj(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len))
    elseif fn == "ϕ"
        return real(ϕ(x; m=m, len=len).*conj(ϕ(x; m=m, len=len)))
    end
end

export ψ
