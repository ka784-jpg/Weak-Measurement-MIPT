using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles
using Random, Distributions

HD = (1/sqrt(2))*[1 1; 1 -1]
ID2 = [1 0; 0 1]
HD = kron(HD,ID2)
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
EG = CNOT*HD

wid = 1
center = 4
Gup = Normal(center,wid/sqrt(2))
Gdown = Normal(-center,wid/sqrt(2))

function Gauss(x,w) 
    G = exp((-x^2)/(2*w^2))/((pi)^(1/4)*(w)^(1/2))
end

function f(x)
        if x == 1
            return [0,1]
        end
        if x == 0 
            return [1,0]
        end
end

# Turns the binary list to the base ten number
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

function splitarb(state,left,right)
    statenew = zeros(right-left+1)
    for i in left:right
        statenew[i-left+1] = Int64(state[i])
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
    for j in range(0,dimen-1)
        l = tob10(splitright(ket[j+1],cut))
        m = tob10(splitleft(ket[j+1],cut))
        A[l+1,m+1] = state[j+1] 
    end
    # Compute the nth Reyni Entropy
    evls = svdvals(A)
        for i in 1:length(evls)
            srn += ((norm(evls[i]))^2n)
        end
        srni = (1/(1-n))*log2(srn)
        return(srni)
end

function reynisp(state,nth,l1,l2)
    # Generate the basis kets, and the Matrix A for SVD
    n = nth
    diml = l2-l1+1
    dimr = nspin-diml
    A = zeros(ComplexF64,2^dimr,2^diml)
    srn = 0 
    for j in range(0,dimen-1)
        l = tob10(cat(splitarb(ket[j+1],1,l1-1), splitarb(ket[j+1],l2+1,nspin),dims =(1,1)))
        m = tob10(cat(splitarb(ket[j+1],l1,l2),dims =(1,1)))
        println([l])
        A[l+1,m+1] = state[j+1] 
    end
    # Compute the nth Reyni Entropy
    evls = svdvals(A)
        for i in 1:length(evls)
            srn += ((norm(evls[i]))^2n)
#            srn += ((norm(evls[i]))^2)*log2((norm(evls[i]))^2)
        end
        srni = (1/(1-n))*log2(srn)
        return(srni)
end

function reynidj(state,nth,l1,l2,l3,l4)
    # Generate the basis kets, and the Matrix A for SVD
    n = nth
    diml = l4+l2-l1-l3+2
    dimr = nspin-diml
    A = zeros(ComplexF64,2^dimr,2^diml)
    srn = 0 
    for j in range(0,dimen-1)
        l = tob10(cat(splitarb(ket[j+1],1,l1-1),splitarb(ket[j+1],l2+1,l3-1),splitarb(ket[j+1],l4+1,nspin),dims=(1,1)))
        m = tob10(cat(splitarb(ket[j+1],l1,l2), splitarb(ket[j+1],l3,l4), dims =(1,1)))
        A[l+1,m+1] = state[j+1] 
    end
    # Compute the nth Reyni Entropy
    evls = svdvals(A)
        for i in 1:length(evls)
           srn += ((norm(evls[i]))^2n)
#            srn += ((norm(evls[i]))^2)*log2((norm(evls[i]))^2)
        end
        srni = (1/(1-n))*log2(srn)
        return(srni)
end

function reyni3(state,nth,l1,l2,l3,l4,l5,l6)
    # Generate the basis kets, and the Matrix A for SVD
    n = nth
    diml = l2+l4+l6-l1-l3-l5+3
    dimr = nspin-diml
    A = zeros(ComplexF64,2^dimr,2^diml)
    srn = 0 
    for j in range(0,dimen-1)
        l = tob10(cat(splitarb(ket[j+1],1,l1-1),splitarb(ket[j+1],l2+1,l3-1),splitarb(ket[j+1],l4+1,l5-1),splitarb(ket[j+1],l6+1,nspin),dims =(1,1)))
        m = tob10(cat(splitarb(ket[j+1],l1,l2), splitarb(ket[j+1],l3,l4), splitarb(ket[j+1],l5,l6), dims = (1,1)))
        A[l+1,m+1] = state[j+1] 
    end
    # Compute the nth Reyni Entropy
    evls = svdvals(A)
        for i in 1:length(evls)
            srn += ((norm(evls[i]))^2n)
#            srn += ((norm(evls[i]))^2)*log2((norm(evls[i]))^2)
        end
        srni = (1/(1-n))*log2(srn)
        return(srni)
end

function mutualinformation(state,nth,l1,l2,l3,l4)
    SA = reynisp(state,nth,l1,l2)
    SB = reynisp(state,nth,l3,l4)
    SAB = reynidj(state,nth,l1,l2,l3,l4)
    return(SA+SB-SAB)
end

function svd_full(psi::Vector{ComplexF64},A::Vector{Int64},L::Int64)

    Schmidt_values = Array{Float64}(undef,2^length(A))
    B = [j for j=1:L if !(j in A)]

    psiSVD = Array{ComplexF64}(undef,2^length(A),2^length(B))

    for RA=CartesianIndices(Tuple(0:1 for j=1:length(A)))
        for RB=CartesianIndices(Tuple(0:1 for j=1:length(B)))
            n = dot(Tuple(RA),2 .^(A .- 1) ) + dot(Tuple(RB),2 .^(B .- 1) ) + 1
            ii = dot(Tuple(RA),2 .^(0:(length(A)-1))) + 1
            jj = dot(Tuple(RB),2 .^(0:(length(B)-1))) + 1
            psiSVD[ii,jj] = psi[n]
        end
    end

    return svdvals(psiSVD)

end

function renyi_entropy(psi::Vector{ComplexF64},A::Vector{Int64},L::Integer,ns::Vector{Number})

    S = Array{Float64,1}(undef,length(ns))

    ps = abs2.(svd_full(psi,A,L))
    for i = eachindex(ns)
        n = ns[i]
        S[i] = Renyi_entropy(n,ps)
    end
    return S
end

function Renyi_entropy(n::Number,ps)
    S = 0
    if n==0
        S = log2(sum(ps .> eps(eltype(ps))))
        #S = log2(sum(ps .> 1e-16))
    elseif n==1
        S = -sum([p*log2(p) for p=ps if 0<p])
    elseif n==Inf
        S = -log2(maximum(ps)^2)
    else
        S = 1/(1-n)*log2(sum([p^n for p=ps if 0<p]))
    end
    return S
end

function mutual_information(psi::Vector{ComplexF64}, A::Vector{Int64},B::Vector{Int64},L::Integer,ns::Vector{Number})

    pA = abs2.(svd_full(psi,A,L))
    pB = abs2.(svd_full(psi,B,L))
    pAB = abs2.(svd_full(psi,unique([A;B]),L))

    SA = Array{Float64}(undef,length(ns))
    SB = Array{Float64}(undef,length(ns))
    SAB = Array{Float64}(undef,length(ns))
    for i = eachindex(ns)
        n = ns[i]
        SA[i] = Renyi_entropy(n,pA)
        SB[i] = Renyi_entropy(n,pB)
        SAB[i] = Renyi_entropy(n,pAB)
    end

    MI = SA .+ SB .- SAB
    return MI
end

function top_mutual_information(psi::Vector{ComplexF64},
        A::Vector{Int64},B::Vector{Int64},C::Vector{Int64},L::Integer,ns::Vector{Number})

    pA = abs2.(svd_full(psi,A,L))
    pB = abs2.(svd_full(psi,B,L))
    pC = abs2.(svd_full(psi,C,L))
    pAB = abs2.(svd_full(psi,unique([A;B]),L))
    pAC = abs2.(svd_full(psi,unique([A;C]),L))
    pBC = abs2.(svd_full(psi,unique([B;C]),L))
    pABC = abs2.(svd_full(psi,unique([A;B;C]),L))

    SA = Array{Float64}(undef,length(ns))
    SB = Array{Float64}(undef,length(ns))
    SC = Array{Float64}(undef,length(ns))
    SAB = Array{Float64}(undef,length(ns))
    SAC = Array{Float64}(undef,length(ns))
    SBC = Array{Float64}(undef,length(ns))
    SABC = Array{Float64}(undef,length(ns))

    for i = eachindex(ns)
        n = ns[i]
        SA[i] = Renyi_entropy(n,pA)
        SB[i] = Renyi_entropy(n,pB)
        SC[i] = Renyi_entropy(n,pC)
        SAB[i] = Renyi_entropy(n,pAB)
        SAC[i] = Renyi_entropy(n,pAC)
        SBC[i] = Renyi_entropy(n,pBC)
        SABC[i] = Renyi_entropy(n,pABC)
    end 
    
    MI = SA .+ SB .+ SC .- SAB .- SAC .- SBC .+ SABC
    return MI

end

function TMI(state,nth,l1,l2,l3,l4,l5,l6)
    SA = reynisp(state,nth,l1,l2)
    SB = reynisp(state,nth,l3,l4)
    SC = reynisp(state,nth,l5,l6)
#    SAB = reynidj(state,nth,l1,l2,l3,l4)
    SAB = reynisp(state,nth,l1,l4)
    SAC = reynidj(state,nth,l1,l2,l5,l6)
#    SBC = reynidj(state,nth,l3,l4,l5,l6)
    SBC = reynisp(state,nth,l3,l6)
#    SABC = reyni3(state,nth,l1,l2,l3,l4,l5,l6)
    SABC = reynisp(state,nth,l1,l6)
    return(SA+SB+SC+SAB+SBC+SAC+SABC)
end

function gate!(statein,i)
    L = nspin
    U = rangate()
    psi_temp = zeros(ComplexF64, 4)
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
    statein .= statein/norm(statein)
end
    
function weakmsman!(statein)
# Weak Measurements
    global xo
    xo = 0
    L = nspin 
    for i=2:L-2
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1)+1
                pup += statein[n]*conj(statein[n])
            end
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
        statein .= statein/norm(statein)
    end
end

function engate!(state, lspin, rspin)
    nl = lspin
    nr = rspin
    dimen = length(state)
    for i=1:dimen
#        if splicee(ket[i],nl,nr) != [0; 0]
#            state[i] = 0 
#        end
    end
    state .= state/norm(state)
    A = zeros(ComplexF64,4,2^(nspin-2))
    f = zeros(Int64,4,2^(nspin-2))
    for j in 1:dimen
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
    state .= state/norm(state)
end

function gate_periodic!(state)
    nl = 2
    nr = L-2
    dimen = length(state)
    U = rangate()
    for i=1:dimen
#        if splicee(ket[i],nl,nr) != [0; 0]
#            state[i] = 0 
#        end
    end
    state .= state/norm(state)
    A = zeros(ComplexF64,4,2^(nspin-2))
    f = zeros(Int64,4,2^(nspin-2))
    for j in 1:dimen
        l = tob10(splicee(ket[j],nl,nr))
        m = tob10(splicer(ket[j],nl,nr))
        A[l+1,m+1] = state[j] 
        f[l+1,m+1] = j 
    end
    A = U*A
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]
            state[f[i,j]] = A[i,j]
        end
    end
    state .= state/norm(state)
end

function runy(nspins, iteration, probability, filename) 
    #Generate Initial State
    prob = probability
    global prob
    nspin = nspins + 2
    L = nspin
    global L
    dimen = 2^nspin
    global nspin
    global dimen
    state = zeros(Complex{Float64}, dimen)
    state[rand(1:dimen)] = 1 
    ket = [zeros(Int64,nspin)];
    global ket
    kettemp = zeros(Int64,nspin)
    global kettemp
    # Initialization conditions for Gates
    timestep = Int64(5*nspins+1)
    MI = zeros(timestep)
    iterations = iteration
    vals = 0
    xp = zeros(timestep)
    n1th = 1
    cut = Int64(nspin/2)
    MIvec = Array{Number,1}([1,2,3])
    Nu_1 = zeros(Float64,4*nspins+1,iteration)
    Nu_2 = zeros(Float64,4*nspins+1,iteration)
    Nu_3 = zeros(Float64,4*nspins+1,iteration)
    dimen = 2^nspin
    global dimen
    for j in 1:dimen-1
        temp = Vector(binaryvec(j))
        push!(ket,temp)
    end
    
    
    for k=1:iterations
        
        state = zeros(Complex{Float64}, dimen)
        state[rand(1:dimen)] = 1
        
        # Run the circuit out to a steady state for t = 2L
        for i=1:Int64(nspins)
            for j=2:2:L-2
                gate!(state,Int64(j))
            end
            weakmsman!(state)            
            for j=3:2:L-3
                gate!(state,Int64(j))
            end
            gate_periodic!(state)
            weakmsman!(state)
        end
        
        # Entangle the ancillas with the end spins
        engate!(state, 1, 2)
        engate!(state, Int64((L-2)/2+2), L)
        
        # Apply the gates and measurements to the system from t=2L to t=10L
        for i=Int64(nspins+1):timestep
            # Even Layer
            for j=2:2:L-2
                gate!(state,Int64(j))
            end
            # Record the Data
            temp = mutual_information(state,[1],[L],L,MIvec)
            Nu_1[Int64(i - nspins),k] += temp[1]/(2)
            Nu_2[Int64(i - nspins),k] += temp[2]/(2) 
            Nu_3[Int64(i - nspins),k] += temp[3]/(2)
            weakmsman!(state)
            # Odd Layer
            for j=3:2:L-3 
                gate!(state,Int64(j)) 
            end
            gate_periodic!(state)
            # Record the Data
            temp = mutual_information(state,[1],[L],L,MIvec)
            Nu_1[Int64(i - nspins),k] += temp[1]/(2)
            Nu_2[Int64(i - nspins),k] += temp[2]/(2)
            Nu_3[Int64(i - nspins),k] += temp[3]/(2)
            weakmsman!(state)
        end
    end
    writedlm(string(filename,"Reyni1"),Nu_1, ", ")
    writedlm(string(filename,"Reyni2"),Nu_2, ", ")
    writedlm(string(filename,"Reyni3"),Nu_3, ", ")
end

runy(8, 1, 0.142, "outputs_8")
