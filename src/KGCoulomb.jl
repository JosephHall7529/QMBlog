include("../../Constants/src/Short.jl")

print("What is the mass of the particle you are dealing with? ")
m0 = parse(Float64,readline())

function ε(n,l; Z=10, m0=m0)
    rest = m0*c^2
    return rest/sqrt(1+((Z*α)^2)/(n-l-1/2+sqrt((l+1/2)^2-(Z*α)^2))^2)
end

β(n,l; kw...) = 2*sqrt(m0^2*c^4 - ε(n,l; kw...)^2)/(ħ*c)

function ρ(r,n,l; kw...)
    return β(n,l; kw...)*r >= 0 ? β(n,l; kw...)*r : error("ρ must be +VE!")
end
Z = 1
l = 0
μ(Z) = sqrt((l+1/2)^2 - (Z*α)^2)
λ(n,l; kw...) = -m0*c^2 ≤ ε(n,l; kw...) ≤ m0*c^2 ? (2Z*α*ε(n,l; kw...))/(ħ*c*β(n,l; kw...)) : error("Not a bund state")

ν(Z) = 1/2 + μ(Z)

function f(r,n,l; N=10, a0=1, kw...)
    try Z == Int catch nothing end == nothing ? Z = 10 : Z = Z
    S = Float64[]
    a = (2Z*α*ε(n,l; kw...))/(ħ*c*β(n,l; kw...))
    b = 2*μ(Z) + 1

    push!(S, a0)

    for i in 2:N
        push!(S, S[i-1]*((a+i-1)*ρ(r,n,l; kw...))/(b+i-1))
    end
    return sum(S)
end

@doc raw"""
    R

We can see examples of this function

"QMBlog/examples/Ex1.10.ipynb"
"""
function R(r,n,l; kw...)
    return ρ(r,n,l; kw...)^(ν(Z))*exp(-ρ(r,n,l; kw...)/2)*f(r,n,l; kw...)
end

export R
