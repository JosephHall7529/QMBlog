using LinearAlgebra

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

# moral
we remove the dependence on p if the particle is localized, no longer requiring it as a variable
now there is also no mechanism to get ϕ(p)
"""
function ψ(x::Float64,t::Float64; p::Float64=x, m::Float64=1., free::String="y", localized::String="n",
            len::Int=1000)
    if free == "y"
        E = sqrt(dot(p.*c,p.*c) + (m*c^2)^2)
        if localized == "n"
            ψ = [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:length(p), j in 1:length(x), t in t]
            return ψ
        elseif localized == "y"
            p = range(-abs(100*x), 100*x, length=len)
            ϕ = rand(len)
            N = 1/(2π*ħ)
            ψ = [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:len, j in 1:length(x), t in t]
            return (1/sqrt(2π*ħ))*(1/len)*sum(ϕ[:] .* ψ for i in 1:len, j in 1:length(x), t in t, dims=1)
        end
    end
end
function ψ(x::Array{Float64,1},t::Array{Float64,1}; p::Array{Float64,1}=x, m::Float64=1., free::String="y",
            localized::String="n", len::Int=1000)
    if free == "y"
        E = sqrt(dot(p.*c,p.*c) + (m*c^2)^2)
        if localized == "n"
            ψ = [exp.(im*(dot(p[i],x[j]) .- E*t) / ħ) for i in 1:length(p), j in 1:length(x), t in t]
            return ψ
        elseif localized == "y"
            p = range(-abs(100*minimum(x)), 100*maximum(x), length=len)
            ϕ = rand(len)
            N = 1/(2π*ħ)
            ψ = [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:len, j in 1:length(x), t in t]
            return (1/sqrt(2π*ħ))*(1/len)*sum(ϕ[:] .* ψ, dims=1)
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
"""
function P(x::Float64,t::Float64; p::Float64=x, m::Float64=1., free::String="y", localized::String="n",
        len::Int=1000)
    return real(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len).*
        conj(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len)))
end
function P(x::Array{Float64,1},t::Array{Float64,1}; p::Array{Float64,1}=x, m::Float64=1., free::String="y",
            localized::String="n", len::Int=1000)
    return real(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len).*
        conj(ψ(x,t; p=p, m=m, free=free, localized=localized, len=len)))
end
