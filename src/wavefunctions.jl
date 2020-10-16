using LinearAlgebra

include("../../Constants/src/Short.jl")

@doc raw"""
the free wave function is given as exp(i(p.x - Et)/ħ)
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

# output

    ψ::Array{Complex,3} - x is horizontal
                          p is vertical
                          t is depth

# moral
currently there is no mechanism by which to get p
"""
function ψ(x::Float64,t::Float64; p::Float64=x, m::Float64=1., free::String="y")
    if free == "y"
        E = sqrt(dot(p.*c,p.*c) + (m*c^2)^2)
        return [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:length(p), j in 1:length(x), t in t]
    end
end
function ψ(x::Array{Float64,1},t::Array{Float64,1}; p::Array{Float64,1}=x, m::Float64=1., free::String="y")
    if free == "y"
        E = sqrt(dot(p.*c,p.*c) + (m*c^2)^2)
        return [exp(im*(dot(p[i],x[j]) - E*t) / ħ) for i in 1:length(p), j in 1:length(x), t in t]
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

# output
P::Array{Complex,3} - x is horizontal
                  p is vertical
                  t is depth

# moral
the probability is 1 everywhere with the old wavefunction defn
"""
function P(x::Float64,t::Float64; p::Float64=x, m::Float64=1., free::String="y")
    return real(ψ(x,t; p=p, m=m, free=free).*conj(ψ(x,t; p=p, m=m, free=free)))
end
function P(x::Array{Float64,1},t::Array{Float64,1}; p::Array{Float64,1}=x, m::Float64=1., free::String="y")
    return real(ψ(x,t; p=p, m=m, free=free).*conj(ψ(x,t; p=p, m=m, free=free)))
end 
