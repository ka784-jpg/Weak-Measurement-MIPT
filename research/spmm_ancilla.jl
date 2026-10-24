using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles
using StatsBase
using Random, Distributions
using Bootstrap

HD = (1/sqrt(2))*[1 1; 1 -1]
ID2 = [1 0; 0 1]
HD = kron(HD,ID2)
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
EG = CNOT*HD
str = 0.45

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
# Divides the Vector into Two Parts at left and right
function splicee(state,left,right)
    statevo = zeros(2)
    statevo[1] = state[left]
    statevo[2] = state[right]
    return round.(Int,statevo)
end

function splicer(state,left,right)
    statevo = zeros(length(state)-2)
    jp = 0 
    for i=1:left-1
        jp = jp + 1
        statevo[jp] = state[i]
    end
    for i=left+1:right-1
        jp = jp + 1
        statevo[jp] = state[i]
    end
    for i=right+1:length(state)
        jp = jp + 1
        statevo[jp] = state[i]
    end
    return round.(Int,statevo)
end

px = [0 1;1 0]
py = [0 -1im; 1im 0]
pz = [1 0; 0 -1]

function Haargate()   
    a = zeros(ComplexF64, 4,4)
    lambda = zeros(ComplexF64, 4,4)
    a = randn(ComplexF64,(4,4))
    Q,R = qr(a)
    lambda = Diagonal(R)/(abs.(Diagonal(R)))
    U = Matrix(Q)*lambda
    return(U)
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
    for j in range(0,dim-1)
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

# Turn the state onto Φ+ Bell State
function engate(state, enspin)
    nl = enspin
    nr = nspin
    for i=1:dim 
#        if splicee(ket[i],nl,nr) != [0; 0]
#            state[i] = 0 
#        end
    end
    state = state/norm(state)
    A = zeros(ComplexF64,4,2^(nspin-2))
    f = zeros(Int64,4,2^(nspin-2))
    for j in 1:dim
        l = tob10(splicee(ket[j],nl,nr))
        m = tob10(splicer(ket[j],nl,nr))
        A[l+1,m+1] = state[j] 
        f[l+1,m+1] = j 
    end
    A = EG*A
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]
            state[f[i,j]] = A[i,j]
        end
    end
    state
end

function entangle(state,i)
    statein = copy(state)
    psi_temp = zeros(ComplexF64, 4)
    U = 1/sqrt(2)*[1 0 1 0; 0 1 0 1; 0 1 0 -1; 1 0 -1 0]

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
    statein
end

function gatean(statein,i)
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
        for nM = 0:4:2^(L-1)-1,nR = 0:1 
            ind1 = nM + nR + 1
            ind2 = nM + nR + 2 + 1
            ind3 = nM + nR + 2^(L-1) + 1
            ind4 = nM + nR + 2^(L-1) + 2 + 1
#            println([ind1])
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

function projmsman(statein)
# Projective Measurements
    L = nspin
    for i=2:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
            p2 = rand()
        # Measure up 
            if p2<Real(pup)
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL +1
                    statein[n] = 0 
                end
        # Measure down 
            else 
                for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                    n = nR + 2^(i)*nL + 2^(i-1)+1
                    statein[n] = 0
                end
            end
        end
        statein = statein/norm(statein)
    end
return(statein)
end

function av(vector)
    k = 0 
    for i=1:length(vector)
        k += vector[i]
    end
    k/length(vector)
end

function weakmsman!(statein)
# Weak Measurements
    L = nspin
    for i=2:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR+2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
            pup = Float64(pup)
            pwup = 1/(2*(1+(str)^2))*(1+(str)^2+2*str*(2*pup-1))
            p2 = rand()
            # Makes the appropriate sign for the Krauss operator depending on up or down measurement
            if p2 < pwup
                sgn = 1        
            else 
                sgn = -1
            end
        # Replaces the spin Down entries
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 1
                statein[n] = ((1-sgn*str)*statein[n])/sqrt(2*(1+str^2))
            end
        # Replaces the spin Up entries
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1)+1
                statein[n] = ((1+sgn*str)*statein[n])/sqrt(2*(1+str^2))
            end
        end
        statein .= statein/norm(statein)
    end
end


function runwkanc(nspins, iteration, probability, filename) 
    #Generate Initial State
    prob = probability
    global prob
    nspin = nspins + 1
    dim = 2^nspin
    global nspin
    global dim
    ket = [zeros(Int64,nspin)];
    global ket
    kettemp = zeros(Int64,nspin)
    global kettemp
    # Initialization conditions for Gates
    timestep = 2*nspin
    iterations = iteration
    x = zeros(iterations)
    vals1 = zeros(iterations)
    vals2 = zeros(iterations)
    vals3 = zeros(iterations)
    n1th = 1
    n2th = 2
    n3th = 3
    cut = nspin-1
    for j in 1:dim-1
        temp = Vector(binaryvec(j))
        push!(ket,temp)
    end
    for k=1:iterations
    state = zeros(ComplexF64,2^nspin)
    state[rand(1:2^nspins)] = 1.0
    state = engate(state,nspin-1)
        for i=1:Int64(timestep/2)
            for j=1:(nspin-1)/2
                state = gatean(state,Int64(2*j))
            end
            for j = 2:(nspin-1)/2+1
                if j == 2 
                    state = gatean(state,Int64(nspin)) 
                end
                if j!=2 
                    state = gatean(state,Int64(2*(j-3)+3))
                end
            end
        end
        for i = Int64(timestep/2+1):timestep
            for j=1:(nspin-1)/2
                state = gatean(state,Int64(2*j))
            end
            if i==timestep
                vals1[k] += reyni(state,n1th,cut)/2 
                vals2[k] += reyni(state,n2th,cut)/2
                vals3[k] += reyni(state,n3th,cut)/2
            end
            weakmsman!(state)
            for j = 2:(nspin-1)/2+1
                if j == 2 
                    state = gatean(state,Int64(nspin)) 
                end
                if j!=2 
                    state = gatean(state,Int64(2*(j-3)+3))
                end
            end
            if i==timestep
                vals1[k] += reyni(state,n1th,cut)/2 
                vals2[k] += reyni(state,n2th,cut)/2
                vals3[k] += reyni(state,n3th,cut)/2
            end
            weakmsman!(state)
        end
    end
    s1 = mean(vals1)
#    s2 = mean(vals2)
#    s3 = mean(vals3)
     writedlm(string(filename,"it_reyni1"),vals1, ", ")
     writedlm(string(filename,"it_reyni2"),vals2, ", ")
     writedlm(string(filename,"it_reyni3"),vals3, ", ")
#     writedlm(string(filename,"reyni1"),s1, ", ")
#     writedlm(string(filename,"reyni2"),s2, ", ")
#     writedlm(string(filename,"reyni3"),s3, ", ")
    return(s1)
end

runwkanc(8,1,.1,"8spinsout1")
