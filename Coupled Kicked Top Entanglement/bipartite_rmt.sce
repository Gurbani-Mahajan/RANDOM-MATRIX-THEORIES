//RMT ANALYSIS OF BIPARTITE SYSTEM
clc
clear

//building the bipartite system rdm
j=[10,11]
nnsd_list=list()
evals_ensemble_list=list()
for l=1:2
j1=j(l)
m=linspace(-j1,j1,2*j1+1) //basis
//building angular momentum operators
//ladder operators
J_plus=zeros(2*j1+1,2*j1+1)
J_minus=zeros(2*j1+1,2*j1+1)
for i=1:length(m)
    element1=sqrt(j1*(j1+1)-m(i)*(m(i)+1))
    element2=sqrt(j1*(j1+1)-m(i)*(m(i)-1))
    if i~= length(m)
    J_plus(i,i+1)=element1 //0 everywhere except at <m|J+|m+1>
    end
    if i ~= 1
    J_minus(i,i-1)=element2 //0 everywhere except at <m|J|m-1>
    end
end


//building Jz,Jx,Jy
J_z=diag(m)
J_x=(J_plus + J_minus)/2
J_y=(J_plus - J_minus)/(2*%i)

//building unitary matrices to evolve kicked top model
K=[3,6,10]
a=0.1 //chosen to create chaotic system to compare to LOE
nnsd_list(l)=list()
evals_ensemble_list(l)=list()
epsi2=input('enter epsi=') 
for k1=1:3
k=K(k1)
nnsd_list(l)(k1)=[]
evals_ensemble_list(l)(k1)=[]
U_f= expm(-%i*(%pi/2)*J_y) //free precession
U_k = expm(-%i*(k/2*j1)*(J_z+a)^2) //torsion
U=U_f*U_k //floquet operator for kicked top

//coupled kicked tops to measure entanglement 
//assuming identical
n=2*j1+1
m=2*j1+1
U_product= kron(U,U) 

U_12=expm(-%i*(epsi2/sqrt(j1*j1))*kron(J_z,J_z)) //interaction operator
U_t=U_product*U_12 //final complete floquet operator
//eigenvectors
[psi_m,eigs_m]=spec(U_t)
evals_ensemble=[]
nnds=[]
for y=1:n*m
    psi=psi_m(:,y) //single eigenstate
        psi_reshape=matrix(psi,n,m)
    //building reduced density matrix
    rdm_1 = psi_reshape*psi_reshape'
    evals_rdm=spec(rdm_1)
    evals_ensemble_list(l)(k1)=[evals_ensemble_list(l)(k1);evals_rdm]
    
    
    //unfolding eigenvalues to find NNSD
    len=length(evals_rdm)
    kicked_eigvals=gsort(evals_rdm,"g","i") //sorting
    
    //removing edge values for bulk stats
    cut=((5/100)*len)
    kicked_evals_bulk=kicked_eigvals(cut+1:$-cut) //removing 5% from edges
    N0=length(kicked_evals_bulk) 
    
    //unfolding using Marchenko-Pastur density
    N_m=2*j1+1
    M_m=2*j1+1
    Q_m=N_m/M_m
    z_max=(1+(1/Q_m)+(2/sqrt(Q_m)))/N_m
    z_min=(1+(1/Q_m)-(2/sqrt(Q_m)))/N_m
    for h=1:N0
        z0=kicked_evals_bulk(h)
        z_vec=linspace(z_min,z0,500)
        mp_density_m=[]
        for g=1:500
            z=z_vec(g)
            mp_density_m(g)=(N_m*Q_m)*sqrt((z_max-z)*(z-z_min))*(1/(2*%pi*z))
        end
        mp_density_m(1)=0 //avoiding singularity
        kicked_evals_unfold(h)=inttrap(z_vec,mp_density_m')
    end
    
    /*
    fit=polyfit(kicked_eigvals,N,2) //interpolating 
    kicked_evals_unfold=polyval(fit,kicked_eigvals) //evaluating interpolating polynomial
    */
    
    //computing nnds
    spacings=[]
    for x = 1:(length(kicked_evals_unfold)-1)
        space= kicked_evals_unfold(x+1)- kicked_evals_unfold(x)
        spacings(x)=space
    end
    spacings=spacings/mean(spacings)
    nnsd_list(l)(k1) = [nnsd_list(l)(k1);spacings]
end
end
end

//spectral density plot
for lsd=1:2
    for ksd=1:3
        scf(ksd)
        clf()
        evals_ensemble=evals_ensemble_list(lsd)(ksd)
        len=length(evals_ensemble)
        histplot(100,real(evals_ensemble))
        //fitting Marchenko-Pastur distribution for Wishart ensemble
        j1_p=j(lsd)
        N=2*j1_p+1
        M=2*j1_p+1
        Q=N/M
        l_max=(1+(1/Q)+(2/sqrt(Q)))/N
        l_min=(1+(1/Q)-(2/sqrt(Q)))/N
        lambda=linspace(l_min,l_max,500)
        for i=1:500
            mp_density(i)=(N*Q)*sqrt((l_max-lambda(i))*(lambda(i)-l_min))/(2*%pi*lambda(i))
         end
         plot(lambda,mp_density)
         title("Kick strength (k)="+string(K(ksd))+" j= "+string(j(lsd))+ ' epsi = '+ string(epsi2))
         legend('RDM Spectral Density of Kicked Top','Marchenko-Pastur Distribution')
         xlabel('Eigenvalues')
         ylabel('Density')
     end
end


//plotting NNS histogram
for lns=1:2
    for kns=1:3
        scf(3+kns)
        clf()
        nnsd=nnsd_list(lns)(kns)
        histplot(250,real(nnsd))
        title("NNSD for j= "+string(j(lns))+" and k= "+string(K(kns))+" epsilon = "+string(epsi2))
        xlabel("Nearest Neighbour Spacing")
        ylabel("Frequency")
        //plotting LOE Wigner surmise to visualise fit
        s0=max(real(nnsd))
        x0=linspace(0,s0,500)
        wigner_loe=(%pi*x0.*exp((-%pi*x0.^2)/4))/2
        plot2d(x0,wigner_loe,style=5)
        legend('histogram of NNDS of RMD','LOE Wigner Surmise')
        h=gca()
        h.data_bounds=[0,0;3,4]
     end
end

        
