# Author Mohamed Kamal AbdElrahman
# FD code for 1D Waveguides for the TE and TM modes


using SparseArrays
using LinearAlgebra
using Plots
# 1D Forward Difference Matrix
const CFloat = Complex{Float64}

function createForwardDifferenceMatrix(Nw::Int)
    d = -ones(CFloat,Nw)
    u = ones(CFloat,Nw-1)
    D = spdiagm(0=> d,1=> u)
    return D
end 
# 1D Backword Difference Matrix
function createBackwardDifferenceMatrix(Nw::Int)
    d = ones(CFloat,Nw)
    l = -ones(CFloat,Nw-1)
    D = spdiagm(0=> d,-1=> l)
    return D
end

function D(s::String,dx ,Lx)
    Nx = trunc(Int,Lx/dx); 
    if s == "f"
    dxf = createForwardDifferenceMatrix(Nx)/dx
    return dxf
    end
    if s == "b"
    dxb = createBackwardDifferenceMatrix(Nx)/dx
    return dxb
    end
end


function getMode(ϵr,λ0,pol::String,Lx,m::Int = 1)
    ϵ0 =  8.85418782e-12
    μ0 = 1.25663706e-6 
    C0 = sqrt(1/ϵ0/μ0)
    k0 = 2pi/λ0
    ω =  k0*C0
    
    N = length(ϵr)
    dx = Lx/N;
    
    
    Df = D("f",dx,Lx)
    Db = D("b",dx,Lx)
    
    T_eps = diagm(0=> ϵ0.*ϵr)
    T_epsxinv = Diagonal((ϵ0 .* ϵr).^-1);
    
    if pol == "TE"
        A = ω^2 .*μ0.*T_eps + Df*Db
    end
    
    if pol == "TM"
        A = ω^2*μ0*T_eps + T_eps*Df*T_epsxinv*Db;
    end
   (vals,vecs) =  eigen(A)
    
    β_m = (vals[end+1-m])^.5;
    
    neff_m  =  β_m/k0

    return(real(neff_m), real(vecs[:,end+1-m]))
end

#Buidling the waveguide, assymetric step index waveguide in this case
nf = 1.5; Lf = 5;
nc = 1.4; Lc = 20;
ns = 1.45;Ls = 20;

Lx = Lf+ Lc + Ls #Computational domain length
Nx = 100; #totalNumberofCells

function generateThreeLayerStepIndexWG(ns, nf,nc,Ls,Lf,Lc,Nx)
    ϵr = ones(Nx);
    dx = (Ls+Lf+Lc)/Nx
    for i = 1:Nx
    if(0<=i*dx<= Ls)
        ϵr[i] = (ns)^2
    elseif (Ls<i*dx<= Lf+Ls)
        ϵr[i] = (nf)^2
    elseif  (Lf+Ls<i*dx<= Lc+Lf+Ls)
         ϵr[i] = (nc)^2
    end
    end
    return ϵr
end

 
ϵr = generateThreeLayerStepIndexWG(ns, nf,nc,Ls,Lf,Lc,Nx)
#########################################################################
#Mode parameters
λ0 = 1.0; 
modeNumber = 1
polarization = "TE";
#########################################################################
#Solving the eigen value problem and getting the mode
mode = getMode(ϵr,λ0,polarization,Lx,modeNumber)
neff = mode[1]
modeProfile  = mode[2]
#########################################################################
#Plotting the field profile
plot(modeProfile, 
linestyle=(:solid), linealpha=1, linewidth=4,title =  "neff = "*repr(neff), label = "")






psi = zeros(10)
psi_0 = ones(10)
objectiveFunction(psi,psi_0)

using Nabla

function calculateModeSensitivty(ϵr,λ0,pol::String,Lx,m::Int = 1)
    # solve the forward problem
    ϵ0 =  8.85418782e-12
    μ0 = 1.25663706e-6 
    C0 = sqrt(1/ϵ0/μ0)
    k0 = 2pi/λ0
    ω =  k0*C0
    
    N = length(ϵr)
    dx = Lx/N;
    
    
    Df = D("f",dx,Lx)
    Db = D("b",dx,Lx)
    
    T_eps = diagm(0=> ϵ0.*ϵr)
    T_epsxinv = Diagonal((ϵ0 .* ϵr).^-1);
    
    if pol == "TE"
        A = ω^2 .*μ0.*T_eps + Df*Db
    end
    
    if pol == "TM"
        A = ω^2*μ0*T_eps + T_eps*Df*T_epsxinv*Db;
    end
   (vals,vecs) =  eigen(A)
    
    lambda_m = (vals[end+1-m]);
    psi_m =  real(vecs[:,end+1-m])

    objectiveFunction(psi,psi_0) =  sum(abs.(psi .- psi_0).^2)

    psi_0 = generatepsi_0(Lx,N) 
    dgdpsi =  ∇(psi_m->objectiveFunction(psi_m, psi_0))(psi_m)
    A_adj =  A- lambda_m*sparse(1.0I,N,N)
    src_adj =(-sparse(1.0I,N,N) + psi_m * psi_m') * dgdpsi[1]  
    psi_m_adj =  A_adj\src_adj
    
    return(psi_m_adj .*  psi_m)
end

calculateModeSensitivty(ϵr,λ0,polarization,Lx,modeNumber)
i = 1
while i<= 1000
    i = i+1;
   ϵr =  ϵr + real.(calculateModeSensitivty(ϵr,λ0,polarization,Lx,modeNumber))
#########################################################################
#Plotting the field profile
mode = getMode(ϵr,λ0,polarization,Lx,modeNumber)
modeProfile  = mode[2]
c = plot( [modeProfile,ϵr ], linestyle=(:solid), linealpha=1, linewidth=4,title =  "neff = "*repr(neff), label = "")
    display(c)
end

plot(ϵr)


function generatepsi_0(L,Nx)
    psi_0 = ones(Nx);
    dx = L/Nx
    for i = 1:Nx
    if(0<=i*dx<= L)
        psi_0[i]= 1-abs(i*dx-L/2)
  
    else  (L/2<i*dx<= L)
         psi_0[i] =  2-2*(1-0)/L*(i*dx)
    end
    end
    return psi_0
end

plot(generatepsi_0(1.0,100))

  psi_0
