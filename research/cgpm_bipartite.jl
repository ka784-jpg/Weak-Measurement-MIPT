using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles
using StatsBase
using Random, Distributions
using Bootstrap
using Statistics

HD = (1/sqrt(2))*[1 1; 1 -1]
ID2 = [1 0; 0 1]
HD = kron(HD,ID2)
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
EG = CNOT*HD
wid = 1
center = 1
Gup = Normal(center,wid/sqrt(2))
Gdown = Normal(-center,wid/sqrt(2))

function Gauss(x,w) 
    G = exp((-x^2)/(2*w^2))/((pi)^(1/4)*(w)^(1/2))
end

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
        elseif n!=Inf 
            srni = (1/(1-n))*log2(srn)
        elseif n==Inf
            srni = -log2(maximum(evls))
        end
        return(srni)
end

function reyniv(state,nth::Vector{Float64},cut)
    # Generate the basis kets, and the Matrix A for SVD
    diml = cut
    dimr = nspin-cut
    A = zeros(ComplexF64,2^dimr,2^diml)
    for j in range(0,dimen-1)
        l = tob10(splitright(ket[j+1],cut))
        m = tob10(splitleft(ket[j+1],cut))
        A[l+1,m+1] = state[j+1] 
    end

    # Create a temp vector to store the Reyni Entropies 
    temp = zeros(Float64,length(nth))
    
    # Compute the Von-Neumann Entropy and the nth Reyni Entropy
    evls = svdvals(A)  
    for i in range(1,length(nth))
        if nth[i]==1
            temp[i] = -sum(((norm(evls[i]))^2)*log2(norm(evls[i])^2) for i=1:length(evls) if 0 < (norm(evls[i]))^2)
        elseif nth[i]!=Inf 
            srn = 0 
            for k in 1:length(evls)
                if 0 < ((norm(evls[k]))^(2*nth[i])) srn += ((norm(evls[k]))^(2*nth[i])) end
            end
            temp[i] = (1/(1-nth[i]))*log2(srn)
        elseif nth[i]==Inf
            temp[i] = -log2(maximum(evls))
        end
    end
    return(temp)
end

function gate!(statein,i)
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
    statein .= statein/norm(statein)
end
    

function projmsm!(statein)
# Strong Measurements
    L = nspin
    for i=1:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR+2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
            pup = Float64(pup)
            p2 = rand()
            # Generates the position of the pointer using P(x) = pup*G(x-λ) + pdown*G(x+λ)
            if p2 < pup
                # Zeros Out the Spin Down Entries
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL +1
                    statein[n] = 0
                end
            else 
                # Zeros Out the Spin Up Entries
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL + 2^(i-1)+1
                    statein[n] = 0
                end
            end 
        end
        statein .= statein/norm(statein)
    end
end

function weakmsm!(statein)
# Projective Measurements
    global xo
    xo = 0
    L = nspin
    for i=1:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
            pup = Float64(pup)
            p2 = rand()
            if p2<=Real(pup)
                global xo
                xo = rand(Gup,1)[1] 
            end
            if p2>Real(pup)   
                global xo
                xo = rand(Gdown,1)[1]
            end
            # Replace Spin Down Entries
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL +1
                    statein[n] =  ((Gauss(xo+center,wid))*statein[n])
                end
        # Replace Spin Up Entries 
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL + 2^(i-1)+1
                    statein[n] =  ((Gauss(xo-center,wid))*statein[n])
                end
        end
        statein = statein/norm(statein)
    end
    return(statein)
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
    timestep = Int64(4*nspin)
    iterations = iteration
    vals = zeros(iteration)
    vals_07 = zeros(Int64(iteration))
    vals_1 = zeros(Int64(iteration))
    vals_2 = zeros(Int64(iteration))
    vals_3 = zeros(Int64(iteration))
    vals_4 = zeros(Int64(iteration))
    vals_5 = zeros(Int64(iteration))
    vals_inf = zeros(Int64(iteration))
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
                gate!(state,Int64(2*j+1))
            end
            if i == timestep
                temp = reyniv(state,[0.7,1,2,3,4,5,Inf],cut)
                vals_07[k] += temp[1]/2
                vals_1[k] += temp[2]/2
                vals_2[k] += temp[3]/2
                vals_3[k] += temp[4]/2
                vals_4[k] += temp[5]/2
                vals_5[k] += temp[6]/2
                vals_inf[k] += temp[7]/2
            end
            weakmsm!(state)
            for j = 0:nspin/2-1
                if j == 0 
                    gate!(state,Int64(nspin)) 
                end
                if j!=0 
                    gate!(state,Int64(2*j))
                end
            end
            if i == timestep
                temp = reyniv(state,[0.7,1,2,3,4,5,Inf],cut)
                vals_07[k] += temp[1]/2
                vals_1[k] += temp[2]/2
                vals_2[k] += temp[3]/2
                vals_3[k] += temp[4]/2
                vals_4[k] += temp[5]/2
                vals_5[k] += temp[6]/2
                vals_inf[k] += temp[7]/2
            end
            weakmsm!(state)
        end
    end
    
    writedlm(string(filename,"_trajectory","reyni0.7"), vals_07, ", ")
    writedlm(string(filename,"_trajectory","reyni1"), vals_1, ", ")
    writedlm(string(filename,"_trajectory","reyni2"), vals_2, ", ")
    writedlm(string(filename,"_trajectory","reyni3"), vals_3, ", ")
    writedlm(string(filename,"_trajectory","reyni4"), vals_4, ", ")
    writedlm(string(filename,"_trajectory","reyni5"), vals_5, ", ")
    writedlm(string(filename,"_trajectory","reyniInf"), vals_inf, ", ")
    
end

runy(12,2,0.05,"testbipartite.out") 
