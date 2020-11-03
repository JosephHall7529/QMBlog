using Roots, Random

include("../../Constants/src/Short.jl")

@doc """
    domain

splits a finite square well potential into 3 regions

# input
- I::Int=2 - which region of the domain are you looking at

## kwargs
- zi::Float64=-50.0 - z axis starting point
- step::Float64=0.1 - step size of domain
- zf::Float64=50.0 - z axis end point
- a::Float64=20.0 - half the width of the potential well

# output
- z::Array{Float64,1}

# example
    julia> domain()
    201-element Array{Float64,1}:
     -10.0
      -9.9
      -9.8
       ⋮
       9.8
       9.9
      10.0
"""
function domain(I::Int=2; zi::Float64=-50.0, step::Float64=0.1, zf::Float64=50.0, a::Float64=20.0)
    if I == 1
        return [i for i in zi:step:-a/2-step]
    elseif I == 2
        return [i for i in -a/2:step:a/2]
    elseif I == 3
        return [i for i in a/2+step:step:zf]
    end
end

@doc """
    potential

creates the square well over the domain

# input
- V0::Float64=50.0 - which region of the domain are you looking at

## kwargs
- zi::Float64=-50.0 - z axis starting point
- step::Float64=0.1 - step size of domain
- - zf::Float64=50.0 - z axis end point
- a::Float64=20.0 - half the width of the potential well

# output
- z::Array{Float64,1}

# example
    julia> potential(;zi=0.,step=1.,zf=5.,a=2.)
    7-element Array{Float64,1}:
     50.0
     50.0
     50.0
      0.0
      0.0
      0.0
      0.0
"""
function potential(V0::Float64=50.0; a::Float64=20.0, zi::Float64=-50.0, step::Float64=0.1, zf::Float64=50.0)

    Z = Float64[]
    for i in 1:3
        l = length(domain(i; a=a, zi=zi, step=step, zf=zf))
        if i == 1 || i == 3
            for j in 1:l
                push!(Z, 0.0)
            end
        elseif i == 2
            for j in 1:l
                push!(Z, V0)
            end
        end
    end
    return Z
end

@doc """
    ψ

produces the dirac 4-vectors for spin up/ down particle/ antiparticles

# input
- z::Array{Float64, 1}=[1.0] - The z domain with which you would like to evaluate the wavefunction over

## kwargs
- E="" - The energy of the particle being evaluated
- m0::Float64=938.0 - The mass of the particle you would like to evaluate
- state::String="" - the state of the dirac particle:
    * "sup" - spin up particle
    * "sdp" - spin down particle
    * "suap" - spin up antiparticle
    * "sdap" - spin down antiparticle
- rel::String="y" - relativistic or not ("y", "n")
- s::Int=42 - Pseudorandom number generator value
- dir::String="f" - which direction is the particle moving
    * "f" - forwards
    * "b" - backwards
- β="" - the width of the guassian distribution

# output
- ψ::Array{Array{Complex{Float64},2},1} - The  4-D matrix for the dirac particle at every point z in your domain

# example
    julia> ψ()
    1-element Array{Array{Complex{Float64},2},1}:
    [0.672663207537872 + 0.33169897380736924im 0.0 + 0.0im 0.00029202784544425176 + 0.00014400272762291917im 0.0004777886170901285 - 0.00012060963495132284im; 0.0 + 0.0im 0.672663207537872 + 0.33169897380736924im 0.00019519682674329236 + 0.00045246750766584097im -0.00029202784544425176 - 0.00014400272762291917im; 0.00029202784544425176 + 0.00014400272762291917im 0.0004777886170901285 - 0.00012060963495132284im 0.672663207537872 + 0.33169897380736924im 0.0 + 0.0im; 0.00019519682674329236 + 0.00045246750766584097im -0.00029202784544425176 - 0.00014400272762291917im 0.0 + 0.0im 0.672663207537872 + 0.33169897380736924im]

# subtle parts to notice
- the function assumes a Gaussian distributed momentum in the x and y directions. Increasing the value of β reduces the effect of the other two dimensions.
-
"""
function ψ(z::Array{Float64, 1}=[1.0]; E="", m0::Float64=938.0, state::String="", rel::String="y", s::Int=42, dir::String="f", β="")

    if β == ""
        β = 100000
    end
    l = length(z)
    p = []

    if E == ""
        E = (m0*c^2)/8
    end

    for i in 1:l
        if dir == "b"
            k = vcat(exp(-(z[i]/β)^2)*rand(Random.seed!(s),2), -v(m0,E))
            push!(p, [k[1]+im*k[2], k[1]-im*k[2], k[3]])
        else
            k = vcat(exp(-(z[i]/β)^2)*rand(Random.seed!(s),2), v(m0,E))
            push!(p, [k[1]+im*k[2], k[1]-im*k[2], k[3]])
        end
    end

    Ω = Array[]
    for i in 1:l
        ω = ComplexF64[]
        push!(ω, 1 + 0im)
        push!(ω, 0 + 0im)
        if rel == "n"
            [push!(ω, 0) for i in 1:2]
        else
            push!(ω, (p[i][3]*c)/(E+m0*c^2))
            push!(ω, (p[i][1]*c)/(E+m0*c^2))
        end

        push!(ω, 0 + 0im)
        push!(ω, 1 + 0im)
        if rel == "n"
            [push!(ω, 0) for i in 1:2]
        else
            push!(ω, (p[i][2]*c)/(E+m0*c^2))
            push!(ω, (-p[i][3]*c)/(E+m0*c^2))
        end

        if rel == "n"
            [push!(ω, 0 + 0im) for i in 1:2]
        else
            push!(ω, (p[i][3]*c)/(E+m0*c^2))
            push!(ω, (p[i][1]*c)/(E+m0*c^2))
        end
        push!(ω, 1 + 0im)
        push!(ω, 0 + 0im)

        if rel == "n"
            [push!(ω, 0 + 0im) for i in 1:2]
        else
            push!(ω, (p[i][2]*c)/(E+m0*c^2))
            push!(ω, (-p[i][3]*c)/(E+m0*c^2))
        end
        push!(ω, 0 + 0im)
        push!(ω, 1 + 0im)
        ω = reshape(ω, 4, 4)
        push!(Ω, ω)
    end

    if state == "sup"
        return [sqrt(Complex((E+m0*c^2)/(2m0*c^2)))*Ω[i][:,1]*exp(im*p[i][3]*z[i]/ħ) for i in 1:l]
    elseif state == "sdp"
        return [sqrt(Complex((E+m0*c^2)/(2m0*c^2)))*Ω[i][:,2]*exp(im*p[i][3]*z[i]/ħ) for i in 1:l]
    elseif state == "sdap"
        return [sqrt(Complex((E+m0*c^2)/(2m0*c^2)))*Ω[i][:,3]*exp(im*p[i][3]*z[i]/ħ) for i in 1:l]
    elseif state == "suap"
        return [sqrt(Complex((E+m0*c^2)/(2m0*c^2)))*Ω[i][:,4]*exp(im*p[i][3]*z[i]/ħ) for i in 1:l]
    else
        return [sqrt(Complex((E+m0*c^2)/(2m0*c^2)))*Ω[i]*exp(im*p[i][3]*z[i]/ħ) for i in 1:l]
    end
end

@doc """
    γ

This is a useful definition made while calculating the relationship between the constants of the total wavefunction.

# input
- E::Float64 - the energy of the particle
- m0::Float64 - the mass of the particle
- V0::Float64 - the potential well magnitude

# output
- γ::Complex

# example
    julia> γ(200.0,938.0,50.0)
    0.946255439029009 + 0.0im
"""
γ(E::Float64, m0::Float64, V0::Float64) = sqrt(Complex(((E-(m0*c^2))*(E-V0+(m0*c^2)))/((E+(m0*c^2))*(E-V0-(m0*c^2)))))

@doc """
    α, β

These matrices determine the relationship between incoming and scattered waves, and the scattered and outgoing waves respectively.

# input
- E::Float64 - the energy of the particle
- m0::Float64 - the mass of the particle
- V0::Float64 - the potential well magnitude
- a::Float64 - half the width of the potential well

# output
- 2x2 Array{Complex{Float64},2}

# example
    julia> map.(α,200., 938., 50., 20.)
    2×2 Array{Complex{Float64},2}:
       0.851445+0.576754im   0.00747969+0.0273958im
     0.00747969-0.0273958im    0.851445-0.576754im

     julia> map.(β, 200., 938., 50., 20.)
    2×2 Array{Complex{Float64},2}:
       0.805685+0.545756im   -0.0070777+0.0259235im
     -0.0070777-0.0259235im    0.805685-0.545756im
"""
α = Function[]
push!(α, (E, m0, V0, a) -> ((γ(E,m0,V0)+1)/2γ(E,m0,V0))*exp(im*(v(m0,E)-v(m0,E-V0))*a/2ħ))
push!(α, (E, m0, V0, a) -> ((γ(E,m0,V0)-1)/2γ(E,m0,V0))*exp(-im*(v(m0,E)+v(m0,E-V0))*a/2ħ))
push!(α, (E, m0, V0, a) -> ((γ(E,m0,V0)-1)/2γ(E,m0,V0))*exp(im*(v(m0,E)+v(m0,E-V0))*a/2ħ))
push!(α, (E, m0, V0, a) -> ((γ(E,m0,V0)+1)/2γ(E,m0,V0))*exp(im*(v(m0,E-V0)-v(m0,E))*a/2ħ))
α = reshape(α, 2, 2)
β = Function[]
push!(β, (E, m0, V0, a) -> ((γ(E,m0,V0)+1)/2)*exp(im*(v(m0,E)-v(m0,E-V0))*a/2ħ))
push!(β, (E, m0, V0, a) -> ((1-γ(E,m0,V0))/2)*exp(im*(v(m0,E)+v(m0,E-V0))*a/2ħ))
push!(β, (E, m0, V0, a) -> ((1-γ(E,m0,V0))/2)*exp(-im*(v(m0,E)+v(m0,E-V0))*a/2ħ))
push!(β, (E, m0, V0, a) -> ((γ(E,m0,V0)+1)/2)*exp(im*(v(m0,E-V0)-v(m0,E))*a/2ħ))
β = reshape(β, 2, 2)

A(C,Cd; E::Float64=(938.0*c^2)/8, m0::Float64=938.0, V0::Float64=50.0, a::Float64=20.0) = (map.(α, E, m0, V0, a) * map.(β, E, m0, V0, a) * [C, Cd])[1]
Ad(C, Cd; E::Float64=(938.0*c^2)/8, m0::Float64=938.0, V0::Float64=50.0, a::Float64=20.0) = (map.(α, E, m0, V0, a) * map.(β, E, m0, V0, a) * [C, Cd])[2]
B(C,Cd; E::Float64=(938.0*c^2)/8, m0::Float64=938.0, V0::Float64=50.0, a::Float64=20.0) = (map.(β, E, m0, V0, a) * [C, Cd])[1]
Bd(C,Cd; E::Float64=(938.0*c^2)/8, m0::Float64=938.0, V0::Float64=50.0, a::Float64=20.0) = (map.(β, E, m0, V0, a) * [C, Cd])[2]

@doc """
    ψd

The wavefunction ψ over the whole domain

# input
- C - an arbitrarily chosen constant by which all other constants defining the wavefuncton are evaluated w.r.t
- Cd - a constant that can be derived by finding the root of the normalisation condition

## kwargs
- E="" - Energy of the particle
- m0::Float64=938.0 - The mass of the particle
- state::String="sup" - the state of the partile
    * "sup" - spin up particle
    * "sdp" - spin down particle
    * "suap" - spin up antiparticle
    * "sdap" - spin down antiparticle
- rel::Float64="y" - relativistic ("y"/ "n")
- a::Float64=20.0 - half the width of the potential well
- s::Float64=42.0 - pseudorandom generator value
- zi::Float64=-50.0 - z axis starting point
- step::Float64=0.1 - step size of domain
- zf::Float64=50.0 - z axis end point
- β="" - the width of the guassian distribution

# output
- Array{Array{Complex{Float64},1},1}

# example
    julia> ψd(0.2,0.2)
    1001-element Array{Array{Complex{Float64},1},1}:
     [-0.13813229634234925 + 0.0im, 0.0 + 0.0im, 0.0 - 0.0001099906960142059im, -9.69290888388749e-16 - 8.253944535837097e-16im]
     [-0.12638465940596072 + 0.0im, 0.0 + 0.0im, 0.0 - 0.00011262161835027968im, -9.800297957906111e-16 - 8.345391125433953e-16im]
     [-0.11437181710378572 + 0.0im, 0.0 + 0.0im, 0.0 - 0.00011501621566375629im, -9.798578996460599e-16 - 8.343927353041114e-16im]
     ⋮
     [-0.20392618459976558 + 0.0im, 0.0 + 0.0im, 0.0 - 9.552401590658669e-5im, -1.7470972131485612e-15 - 1.4877312547541704e-15im]
     [-0.1936355792485479 + 0.0im, 0.0 + 0.0im, 0.0 - 9.947822036020461e-5im, -1.5015163872000826e-15 - 1.2786082205107243e-15im]
     [-0.18293864930963696 + 0.0im, 0.0 + 0.0im, 0.0 - 0.00010322367987640355im, -1.2837024403800632e-15 - 1.0931299231574157e-15im]
"""
function ψd(C, Cd; V0::Float64=50.0, E="", m0::Float64=938.0, state::String="sup", rel="y", a::Float64=20.0, s::Int=42, zi::Float64=-50.0, step::Float64=0.1, zf::Float64=50.0, β="")

    if E == ""
        E = (m0*c^2)/8
    end
    if β == ""
        β = 1000000
    end

    x = domain(1; a=a, zi=zi, step=step, zf=zf)
    ψd = A(C,Cd; E=E, m0=m0, V0=V0, a=a).*ψ(x; E=E, m0=m0, state=state, rel=rel, s=s, β=β).+
        Ad(C,Cd; E=E, m0=m0, V0=V0, a=a).*ψ(x; E=E, m0=m0, state=state, rel=rel, s=s, dir="b", β=β)

    y = domain(2; a=a, zi=zi, step=step, zf=zf)
    ψd = cat(ψd, B(C,Cd; E=E, m0=m0, V0=V0, a=a).*ψ(y; E=E-V0, m0=m0, state=state, rel=rel, s=s, β=β).+
        Bd(C,Cd; E=E, m0=m0, V0=V0, a=a).*ψ(y; E=E-V0, m0=m0, state=state, rel=rel, s=s, dir="b", β=β), dims=1)

    z = domain(3; a=a, zi=zi, step=step, zf=zf)
    ψd = cat(ψd, C.*ψ(z; E=E, m0=m0, state=state, rel=rel, s=s, β=β).+
        Cd.*ψ(z; E=E, m0=m0, state=state, rel=rel, s=s, dir="b", β=β), dims=1)

    return ψd
end

@doc """
    Const


"""
function Const(C=1e-4; Cd=-0, V0::Float64=50.0, E="", m0::Float64=938.0, state::String="sup", rel="y", a::Float64=20.0, s::Int=42, zi::Float64=-50.0, step::Float64=0.1, zf::Float64=50.0, β="")

    if E == ""
        E = (m0*c^2)/8
    end
    if β == ""
        β = 1000000
    end

    l = length(zi:step:zf)
    an(Cd) = ψd(C, Cd; V0=V0, E=E, m0=m0, state=state, rel=rel, a=a, s=s, zi=zi, step=step, zf=zf, β=β)

    function Φ(Cd)
        Δ = an(Cd)
        I = 0.0im
        for i in 1:l-1
            I += sum(conj(Δ[i]) .* Δ[i])*step
        end
        return Real(I) - 1.
    end

    if Cd == -0

        return find_zero(Φ, (-10.0im,10.0), Order5())
    else
        return Φ(Cd)
    end
end

function dirac(C=1e-4; V0::Float64=50.0, E="", m0::Float64=938.0, state::String="sup", rel::String="y", a::Float64=20.0, s::Int=42, zi::Float64=-50.0, step::Float64=0.1, zf::Float64=50.0, bilinear::String="s", β="")

    if E == ""
        E = (m0*c^2)/8
    end
    if β == ""
        β = 1000000
    end

    println("domain is: ", zi, ":", step, ":", zf)
    Cd = Const(C; V0=V0, E=E, m0=m0, state=state, rel=rel, a=a, s=s, zi=zi, step=step, zf=zf, β=β)

    if bilinear == ""
        return ψd(C, Cd; V0=V0, E=E, m0=m0, state=state, rel=rel, a=a, s=s, zi=zi, step=step, zf=zf, β=β)
    elseif bilinear == "s"
        LL = ψd(C, Cd; V0=V0, E=E, m0=m0, state=state, rel=rel, a=a, s=s, zi=zi, step=step, zf=zf, β=β)
        L = Float64[]
        for (i,j) in enumerate(zi:step:zf)
            push!(L, real(sum(conj(LL[i]) .* LL[i])))
        end
        return L
    end
end

function plotdirac(C=1e-4; V0::Float64=50.0, E="", m0::Float64=938.0, state::String="sup", rel::String="y", a::Float64=20.0, s::Int=42, zi::Float64=-50.0, step::Float64=0.1, zf::Float64=50.0, bilinear::String="s", β="")

    if β == ""
        β = 1000000
    end

    D = dirac(C; V0=V0, E=E, m0=m0, state=state, rel=rel, a=a, s=s, zi=zi, step=step, zf=zf, bilinear=bilinear, β=β)

    upper = maximum(D)
    pot = potential(V0; zi=zi, step=step, zf=zf, a=a).*upper

    P = plot(z, D, background_color="darkgrey", legend=false)
    plot!(z, pot)

    return P
end

export dirac, potential
