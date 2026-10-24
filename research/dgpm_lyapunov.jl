# Ceff from the measurement record
using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles
using StatsBase
using Random, Distributions
using Bootstrap
using SpecialFunctions

HD = (1/sqrt(2))*[1 1; 1 -1]
ID2 = [1 0; 0 1]
HD = kron(HD,ID2)
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
EG = CNOT*HD
wid = 10000
cen = 10000
Gup = Normal(cen,wid/sqrt(2))
Gdown = Normal(-cen,wid/sqrt(2))
no = 1
ep = 10.0^(-no) 

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
            srn += ((norm(evls[i]))^2n)
        end
        srni = (1/(1-n))*log2(srn)
        return(srni)
end


function gate!(statein,i)
    L = nspin
    U = rangate()
    psi_temp = zeros(ComplexF64, 4,2)
    if i!=nspin
        for nL=0:2^(L-1-i)-1, nR=0:2^(i-1)-1
            ind1 = nR + nL*2^(i+1) + 1
            ind2 = nR + nL*2^(i+1) + 2^(i-1) + 1
            ind3 = nR + nL*2^(i+1) + 2^(i) + 1
            ind4 = nR + nL*2^(i+1) + 2^(i) + 2^(i-1) + 1
            psi_temp[1,:] .= statein[ind1,:]
            psi_temp[2,:] .= statein[ind2,:]
            psi_temp[3,:] .= statein[ind3,:]
            psi_temp[4,:] .= statein[ind4,:]
            psi_temp .= U*psi_temp
            statein[ind1,:] .= psi_temp[1,:]
            statein[ind2,:] .= psi_temp[2,:]
            statein[ind3,:] .= psi_temp[3,:]
            statein[ind4,:] .= psi_temp[4,:]
        end
    else 
        for n = 0:(2^(L-2)-1)
            ind1 = 2n + 1
            ind2 = 2n + 1 + 1
            ind3 = 2n + 2^(L-1) + 1
            ind4 = 2n + 2^(L-1) + 1 + 1
            psi_temp[1,:] .= statein[ind1,:]
            psi_temp[2,:] .= statein[ind2,:]
            psi_temp[3,:] .= statein[ind3,:]
            psi_temp[4,:] .= statein[ind4,:]
            psi_temp .= U*psi_temp
            statein[ind1,:] .= psi_temp[1,:]
            statein[ind2,:] .= psi_temp[2,:]
            statein[ind3,:] .= psi_temp[3,:]
            statein[ind4,:] .= psi_temp[4,:]
        end
    end
end

function gsproj!(state,num_vecs,fe)
	q, r = qr(state)
	state .= Matrix(q)
	diagr = abs.(diag(r)).^2
	for i in 1:num_vecs
	    if diagr[i] < eps(Float64)
	        diagr[i] = 0.0
                state[:,i] .= 0.0
	    end
	end
    fe[t:end,1] .+= -log.(diagr)[1]
    fe[t:end,2] .+= -log.(diagr)[2]
end


function weakmsm!(statein,fe)
# Weak Measurements
    global xo
    xo = 0
    L = nspin
    for i=1:L
        p1 = rand()
        pup = Float64(0) 
        if p1<prob
            # Go through all the 1's at the i slot, and sum the up amplitudes
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1) + 1
                pup += statein[n,1]*conj(statein[n,1])
            end
            p2 = rand()
            if p2<=Real(pup)
                global xo
                xo = round(rand(Gup,1)[1],digits = no) 
            end
            if p2>Real(pup)   
                global xo
                xo = round(rand(Gdown,1)[1],digits = no)
            end
            # Replace Spin Down Entries
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 1
                statein[n,:] .= sqrt((erf((2*(xo+cen)+ep)/(2*wid))-erf((2*(xo+cen)-ep)/(2*wid)))/(2))*statein[n,:]   
            end
            # Replace Spin Up Entries 
            for nL=0:(2^(L-i)-1), nR=0:(2^(i-1)-1)
                n = nR + 2^(i)*nL + 2^(i-1)+1
                statein[n,:] .= sqrt((erf((2*(xo-cen)+ep)/(2*wid))-erf((2*(xo-cen)-ep)/(2*wid)))/(2))*statein[n,:]
            end
            fe[t:end,1] .+= -log(norm(statein[1:end,1])^2)
            fe[t:end,2] .+= -log(norm(statein[1:end,2])^2)
        end
        statein[1:end,1] .= statein[1:end,1]/norm(statein[1:end,1])
        statein[1:end,2] .= statein[1:end,2]/norm(statein[1:end,2])
    end
end

function runy(nspins, probability, iterations, filename) 
    global state
    #Generate Initial State
    prob = probability
    global prob
    nspin = nspins
    dimen = 2^nspin
    global nspin
    global dimen
    ket = [zeros(Int64,nspin)];
    global ket
    kettemp = zeros(Int64,nspin)
    global kettemp
    # Initialization conditions for Gates
    timestep = Int64(2*nspin + 12*nspin)
    global timestep
    vals1 = zeros(timestep)
    x = zeros(2*timestep)
    cut = Int64(nspin/2)
    for j in 1:dimen-1
        temp = Vector(binaryvec(j))
        push!(ket,temp)
    end
    for k=1:iterations
        FET = zeros(ComplexF64,2*timestep,2)
        global FET
        FET_write = zeros(ComplexF64,2*timestep,3)
        global FET_write
        t = 0
        global t
        # Generate a new random product state for each trajectory
        state = zeros(Complex{Float64}, dimen, 2)
        state[rand(1:dimen),1] = 1
        state .= Matrix(qr(state).Q)
        for i=1:timestep
            for j=0:nspin/2-1
               gate!(state,Int64(2*j+1))
            end
            t = t + 1
            x[t] = t
            weakmsm!(state,FET)
            gsproj!(state,2,FET)
            for j = 0:nspin/2-1
                if j == 0 
                    gate!(state,Int64(nspin)) 
                end
                if j!=0 
                    gate!(state,Int64(2*j))
                end
            end
            t = t + 1
            x[t] = t
            weakmsm!(state,FET)
            gsproj!(state,2,FET)
        end
        FET_write[1:end,1] = FET[:,1]
        FET_write[1:end,2] = FET[:,2] .- FET[:,1]
        FET_write[1:end,3] = (FET[:,2] .- FET[:,1]).^2
        writedlm(string(filename,"iterat=",k,"_FET"), FET_write/iterations, ", ") 
    end
end

runy(8,0.194,1,"outputs")

