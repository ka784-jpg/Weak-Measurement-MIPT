using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles
using StatsBase
using Random, Distributions
#using Bootstrap
using Statistics

function Gauss(x,w) 
    G = exp((-x^2)/(2*w^2))/((pi)^(1/4)*(w)^(1/2))
end
wid = 1
mean = 1
Gup = Normal(mean,wid/sqrt(2))
Gdown = Normal(-mean,wid/sqrt(2))

function twogaussian(x,mean,w,prob1,prob2)
    F = prob1*Gauss(x-mean,w) + prob2*Gauss(x+mean,w)
end

function f(x)
        if x == 1
            return [0,1]
        end
        if x == 0 
            return [1,0]
        end
end

# Turns the binary list to the base ten number number
function tob10(vecs)
    join(vecs)
    parse(Int, join(vecs); base = 2)
end

# Splits the Vector
function splitleft(state, dimstate)
    statenew = zeros(dimstate)
    for i in 1:dimstate
        statenew[i] = Int64(state[i])
    end
    return round.(Int,statenew)
end

function splitright(state,dimstate)
    statenew = zeros(length(state)-dimstate)
    for i in dimstate+1:length(state)
        statenew[i-dimstate] = Int64(state[i])
    end
    return round.(Int,statenew)
end

function ranm()
    alp = rand(Float64)*2*pi
    psi = rand(Float64)*2*pi
    chi = rand(Float64)*2*pi
    eps = rand(Float64)
    phi = asin(sqrt(eps))
    return exp(1im*alp)*[exp(1im*psi)*cos(phi) exp(1im*chi)*sin(phi); -1*exp(-1im*chi)*sin(phi) exp(-1im*psi)*cos(phi)]
end

function rangate()
    phig = rand(Float64)*2*pi
    J = rand(Float64)*2*pi
    up = ranm()
    um = ranm()
    vp = ranm()
    vm = ranm()
    return exp(1im*phig)*kron(up,um)*exp(-1im*(pi/4)*(kron(px,px)+kron(py,py))-1im*J*kron(pz,pz))*kron(vp,vm)
end

function splitarb(state,left,right)
    statenew = zeros(right-left+1)
    for i in left:right
        statenew[i-left+1] = Int64(state[i])
    end
    return round.(Int,statenew)
end

px = [0 1;1 0]
py = [0 -1im; 1im 0]
pz = [1 0; 0 -1]

# Turns the binary list to the base ten number number
function tob10(vecs)
    join(vecs)
    parse(Int, join(vecs); base = 2)
end

function binaryvec(n)
    ket0 = string(n, base=2, pad = nspin)
    for k in range(1,length(ket0))      
        kettemp[k] = parse(Int64,ket0[k:k]);
    end
    return kettemp
end

function reyni(state,nth,cut)
    # Generate the basis kets, and the Matrix A for SVD
    n = nth
    diml = cut
    dimr = nspin-cut
    A = zeros(ComplexF64,2^dimr,2^diml)
    srn = 0 
    for j in range(0,dimen-1)
        l = tob10(splitright(ket[j+1],cut))
        m = tob10(splitleft(ket[j+1],cut))
        A[l+1,m+1] = state[j+1] 
    end

# Compute the Von-Neumann Entropy and the nth Reyni Entropy
    evls = svdvals(A)  
        for i in 1:length(evls)
            if 0 < ((norm(evls[i]))^2n) srn += ((norm(evls[i]))^2n) end
        end
        if n==1
            srni = -sum(((norm(evls[i]))^2)*log2(norm(evls[i])^2) for i=1:length(evls) if 0 < (norm(evls[i]))^2)
        else
            srni = (1/(1-n))*log2(srn)
        end
        return(srni)
end

function gate(statein,i)
    L = nspin
    U = rangate()
    psi_temp = zeros(ComplexF64, 4)
    if i!=nspin
        for nL=0:2^(L-1-i)-1, nR=0:2^(i-1)-1
            ind1 = nR + nL*2^(i+1) + 1
            ind2 = nR + nL*2^(i+1) + 2^(i-1) + 1
            ind3 = nR + nL*2^(i+1) + 2^(i) + 1
            ind4 = nR + nL*2^(i+1) + 2^(i) + 2^(i-1) + 1
            psi_temp[1] = statein[ind1]
            psi_temp[2] = statein[ind2]
            psi_temp[3] = statein[ind3]
            psi_temp[4] = statein[ind4]
            psi_temp = U*psi_temp
            statein[ind1] = psi_temp[1]
            statein[ind2] = psi_temp[2]
            statein[ind3] = psi_temp[3]
            statein[ind4] = psi_temp[4]
        end
    else 
        for n = 0:(2^(L-2)-1)
            ind1 = 2n + 1
            ind2 = 2n + 1 + 1
            ind3 = 2n + 2^(L-1) + 1
            ind4 = 2n + 2^(L-1) + 1 + 1
            psi_temp[1] = statein[ind1]
            psi_temp[2] = statein[ind2]
            psi_temp[3] = statein[ind3]
            psi_temp[4] = statein[ind4]
            psi_temp = U*psi_temp
            statein[ind1] = psi_temp[1]
            statein[ind2] = psi_temp[2]
            statein[ind3] = psi_temp[3]
            statein[ind4] = psi_temp[4]
        end
    end
    return(statein)
end
    

function weakmsm(statein)
statec = copy(statein)
# Weak Measurements
    L = nspin
    for i=1:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR+2^(i)*nL + 2^(i-1)+1
                pup += statec[n]*conj(statec[n])
            end
            pup = Float64(pup)
            p2 = rand()
# Generates the position of the pointer using P(x) = pup*G(x-λ) + pdown*G(x+λ)
            if p2 < pup
                xo = rand(Gup,1)[1]
            else 
                xo = rand(Gdown,1)[1]
            end
            
        # Replaces the spin Down entries
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL +1
                statec[n] = (Gauss(xo+mean,wid))*statec[n]
            end
        # Replaces the spin Up entries
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1)+1
                statec[n] = (Gauss(xo-mean,wid))*statec[n]
            end
        end
        statec = statec/norm(statec)
    end
    return(statec)
end

function runy(nspins, iteration, probability, filename) 
    #Generate Initial State
    nspin = nspins
    dimen = 2^nspin
    prob = probability
    global nspin
    global dimen
    global prob
    state = zeros(Complex{Float64}, dimen)
    qi = Int64(rand(1:dimen))
    state[qi] = 1
    ket = [zeros(Int64,nspin)];
    global ket
    kettemp = zeros(Int64,nspin)
    global kettemp
    # Initialization conditions for Gates
    timestep = Int64(100)
    iterations = iteration
    vals = zeros(iteration)
    vals1 = zeros(Int64(timestep-nspin/2+1))
    n1th = 1
    for j in 1:dimen-1
        temp = Vector(binaryvec(j))
        push!(ket,temp)
    end
    for k=1:iterations
        cut = Int64(nspin/2)
        sr = rand(1:dimen)
        state = zeros(Complex{Float64}, dimen)
        state[sr] = 1
        for i=1:timestep
            for j=0:nspin/2-1
                gate(state,Int64(2*j+1))
            end
            state = weakmsm(state)
            for j = 0:nspin/2-1
                if j == 0 
                    gate(state,Int64(nspin)) 
                end
                if j!=0 
                    gate(state,Int64(2*j))
                end
            end
            state = weakmsm(state)
            if i > nspin/2 - 1 
                vals1[Int64(i-nspin/2+1)] = reyni(state,n1th,cut)
            end
        end
        vals[k] = var(vals1)
    end
    writedlm(filename, Statistics.mean(vals), ", ")
    writedlm(string(filename,"trajectory"), vals,", ")
end

runy(12,250,0.02,"/scratch/ka784/12spin_") 
