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



plot(generatepsi_0(1.0,100))

  psi_0

using LinearAlgebra
function formTE_EigenProblem(ϵr,λ0,Lx)
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
    A = ω^2 .*μ0.*T_eps + Df*Db
    return(A)
end
function Calculate_α(oldRelativePermitivity,objFuncGrad,∇ϵr_max = 0.1,γ = -0)
    max_dJdϵ = findmax(real.(objFuncGrad))[1]
    max_oldRelativePermitivity = findmax(real.(oldRelativePermitivity))[1]
    α = 1/max_dJdϵ * (∇ϵr_max - γ * max_oldRelativePermitivity)
    return α
end

Nx = 100;
ϵr = generateThreeLayerStepIndexWG(1, nf,1,Ls,Lf,Lc,Nx)
objectiveFunction(psi,psi_0) = sum(abs.(psi .- psi_0).^2)
function generatepsi_0(L,Nx)
    psi_0 = ones(Nx);
    dx = L/Nx
    for i = 1:Nx
        psi_0[i]=  exp(-((i*dx - .45L)/.05L)^2)+ exp(-((i*dx - .55L)/.05L)^2)
    end
    return psi_0
end
a = ((generatepsi_0(Lx,Nx)' * generatepsi_0(Lx,Nx)))^.5
x_0 = (1/a)*generatepsi_0(Lx,Nx)
m = 1


ϵ0 =  8.85418782e-12
μ0 = 1.25663706e-6 
C0 = sqrt(1/ϵ0/μ0)
k0 = 2pi/λ0
ω =  k0*C0

i = 1
for i   = 1:1000
A = formTE_EigenProblem(ϵr,λ0,Lx);
(vals,vecs) =  eigen(A);
x_m = real.(vecs[:,end+1-m])
λ_m = (vals[end+1-m]);

A_adj = A' - λ_m * sparse(1.0I,Nx,Nx);
dgdx =  ∇(x_m->objectiveFunction(x_m, x_0))(x_m)[1]
src_adj = x_m' * dgdx * x_m - dgdx;
beta_m = A_adj\src_adj
sensitivity = real.((ω^2 .*μ0.* (beta_m ) .*  x_m))
step = Calculate_α(ϵr,sensitivity,.0001)
ϵr[25:75] = ϵr[25:75] .-  step .* sensitivity[25:75]
   for i=1:Nx
        if ϵr[i] >= 2
            ϵr[i] = 2
        end
         if ϵr[i] <= 1.95
            ϵr[i] = 1.95
        end
            
    end
end

p0 = plot([x_m,x_0])
p1 = plot([ϵr])
plot(p0,p1)
#plot(sensitivity)

plot(real.(vecs[:,end-1:end+1-m]))


