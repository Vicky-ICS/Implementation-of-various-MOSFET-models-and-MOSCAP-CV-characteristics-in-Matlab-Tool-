q = 1.6e-19;
eps_0= 8.85e-12;
vt=26e-3;

Nsub= -3e17*1e6; %-ve(nmos) +ve(pmos)
k_si=12;
ni=1.5e10*1e6;
Eg=1.1*q;
Lsub= 250e-9;
eps_si= k_si*eps_0;
chi_si=4.05*q;

if Nsub<0
    psub= abs(Nsub); nsub= ni^2/abs(Nsub);
elseif Nsub>0
    psub= ni^2/Nsub; nsub = Nsub;
else
    psub = ni; nsub= ni;
end

%oxide
tox=2.2*1e-9;
k_ox=4;
Eg_ox=9*q;
eps_ox=k_ox*eps_0;
Ec_off= 3.1*q;

%metal related
tm=10e-9;
phi_m= chi_si/q;%Eg/q
phi_b= -sign(Nsub)*vt*log(abs(Nsub)/ni); 
phi_s=chi_si/q+ Eg/(2*q)+ phi_b;
Vfb=phi_m-phi_s;

%solver
dx= 0.1e-9;
NRmax=100;
tol=1e-6;

%mesh 
xox=-tox:dx:0;
xsi= dx:dx:Lsub;
x=[xox,xsi];
N= length (x);
iox= x<=0;
isi=x>0;
i0=find (x==0);

%Newton  Raphson Matrices
ki= k_ox*iox+ k_si*isi;
kip= 0.5*(ki + circshift(ki,-1));
kim= 0.5*(ki + circshift(ki,1));
kipm=kip+kim;
A= diag(kipm,0)- diag (kip(1:end-1),1)- diag (kim(2:end),-1);

 A(1,:)=0; A(1,1)=1;
 b=zeros(N,1);
 
 psi= zeros(N,1);
 
 %fixed oxide charges
 
 rho_ox= zeros(size(xox));
 %rho_ox(round (length(xox)/2))=0*1e12*1e4/dx*q;
 rho_ox(end)= 2e12*1e4/dx*q;
 rho_oxfx= rho_ox;
 
 
 %Interface Charges
 %dite=1e14*1e4/q;
 %ditm= 1e12*1e4/q;
 %E0=Eg/(2*log(2*dite/ditm));
 %ditf=@(E) 1.*(dite.*exp(-E/E0)+ dite*exp(-(Eg-E)/E0));
 %input
 ditf=(0*1e12*1e4/Eg);
 
 
 Vgs= (-2.5:50e-3:2.5);
 Vwiggle=10e-3;
 Vright=0;
 sweep_type=1;% 0= LFCV,1=HFCV,2=Fast Sweep
 dit_dpsi=0.1e-3;
 
 Caps=[];
 psi_ss=[];
 
 %Voltage Sweep
 
     for Vg= Vgs
      
     
     
     V1= Vg-Vfb;
     V2=V1+Vwiggle;
     
     Qsub=[];
     
     Vleft=V1;
     for i=1:NRmax
          %figure(1); plot(x,psi);hold on;
          %Generate b
          
          psi_si= psi(isi);
          psi_s= psi(i0);
          E1= q*psi_s+Eg/2-q*phi_b;
          E2=Eg/2;
          rho_dit= (q/dx)*ditf*(E2-E1);
          rho_dit_dpsi=(q/dx)*ditf*(E2-E1-(q*dit_dpsi));
          delrho_dit = 1/dit_dpsi* (rho_dit_dpsi-rho_dit);
         
                   
          p= psub*exp(-psi_si/vt);
          n= nsub*exp(psi_si/vt);
          delp= -1/vt*p;
          deln=1/vt*n;
          
          if sweep_type==2
              minc=0;
              delminc=0;
              ditq=0;
              ditdq=0;
              if sign(Nsub)==1
                  majc=n;
                  delmajc=deln;
              elseif sign (Nsub)==-1
                  majc=p;
                  delmajc=delp;
              end
          else
              ditq=rho_dit;
              ditdq=delrho_dit;
              if sign (Nsub)==1
                  majc=n;
                  delmajc=deln;
                  minc=p;
                  delminc=delp;
              elseif sign(Nsub)==-1
                  majc=p;
                  delmajc=delp;
                  minc=n;
                  delminc=deln;
              end
          end
          rho_ox(end)= rho_oxfx(end)+ditq;
          rho_si= q*(Nsub+sign (Nsub)*(minc-majc));
          
          rho= [rho_ox';rho_si];
          b=rho*dx^2/eps_0;
          b(1)= Vleft;
          b(N)= Vright;
          f= A*psi-b;
          
          %Jacobian Calculation
          delrho_ox= zeros(size(rho_ox));
          delrho_ox(end)= delrho_ox(end)+ ditdq;
          
          delrho_si= q*(sign(Nsub)*(delminc-delmajc));
          
          delrho= [delrho_ox';delrho_si];
          delb=delrho*dx^2/eps_0;
          delb(1)=0;
          delb(N)=0;
          J= A-diag(delb);
          dV= -J\f;
          
          if max(abs(dV))<tol
              break;
          end
          psi=psi+dV;
     end
     
     Q1= sum(rho_si)*dx+ ditq*dx;
  
     minc1=minc;
     ditq1=ditq;
     
     %Wiggle point
     Vleft = V2;
     for i=1:NRmax
         
         %figure(2);plot(x,psi);hold on;
         
         %generate b
         psi_si= psi(isi);
          psi_s= psi(i0);
          E1= q*psi_s+Eg/2-q*phi_b;
          E2=Eg/2;
          rho_dit= (q/dx)*ditf*(E2-E1);
          rho_dit_dpsi=(q/dx)*ditf*(E2-E1-(q*dit_dpsi));
          delrho_dit = 1/dit_dpsi* (rho_dit_dpsi-rho_dit);
          p= psub*exp(-psi_si/vt);
          n= nsub*exp(psi_si/vt);
          delp= -1/vt*p;
          deln=1/vt*n;
          if sweep_type==2
              minc=0;
              delminc=0;
              ditq=0;
              ditdq=0;
              if sign(Nsub)==1
                  majc=n;
                  delmajc=deln;
              elseif sign (Nsub)==-1
                  majc=p;
                  delmajc=delp;
              end
          
          elseif sweep_type==1
              minc=minc1;
              delminc=0;
              ditq=ditq1;
              ditdq=0;
              if sign(Nsub)==1
                  majc=n;
                  delmajc=deln;
              elseif sign (Nsub)==-1
                  majc=p;
                  delmajc=delp;
              end
           elseif sweep_type==0
               ditq=rho_dit;
               ditdq=delrho_dit;
              if sign (Nsub)==1
                  majc=n;
                  delmajc=deln;
                  minc=p;
                  delminc=delp;
              elseif sign(Nsub)==-1
                  majc=p;
                  delmajc=delp;
                  minc=n;
                  delminc=deln;
              end
          end
          rho_ox(end)= rho_oxfx(end)+ditq;
          rho_si= q*(Nsub+sign (Nsub)*(minc-majc));
          
          rho= [rho_ox';rho_si];
          b=rho*dx^2/eps_0;
          b(1)= Vleft;
          b(N)= Vright;
          f= A*psi-b;
          
          %Jacobian Calculation
          delrho_ox= zeros(size(rho_ox));
          delrho_ox(end)= delrho_ox(end) + ditdq;
          delrho_si= q*(sign(Nsub)*(delminc-delmajc));
          
          delrho= [delrho_ox';delrho_si];
          delb=delrho*dx^2/eps_0;
          delb(1)=0;
          delb(N)=0;
          J= A-diag(delb);
          dV= -J\f;
          
          if max(abs(dV))<tol
              break;
          end
          psi=psi+dV;
     end
     
     Q2= sum(rho_si)*dx+ditq*dx;
     
     psi_ss=[psi_ss psi_s];
     Cap= 1/Vwiggle*(Q1-Q2);
     Caps=[Caps Cap];
     end
      
  
 
    figure(4);
    plot(Vgs,Caps*10^2)
    hold on
    
    q = 1.6e-19;
eps_0= 8.85e-12;
vt=26e-3;

Nsub= -3e17*1e6; %-ve(nmos) +ve(pmos)
k_si=12;
ni=1.5e10*1e6;
Eg=1.1*q;
Lsub= 250e-9;
eps_si= k_si*eps_0;
chi_si=4.05*q;

if Nsub<0
    psub= abs(Nsub); nsub= ni^2/abs(Nsub);
elseif Nsub>0
    psub= ni^2/Nsub; nsub = Nsub;
else
    psub = ni; nsub= ni;
end

%oxide
tox=2.2*1e-9;
k_ox=4;
Eg_ox=9*q;
eps_ox=k_ox*eps_0;
Ec_off= 3.1*q;

%metal related
tm=10e-9;
phi_m= chi_si/q;%Eg/q
phi_b= -sign(Nsub)*vt*log(abs(Nsub)/ni); 
phi_s=chi_si/q+ Eg/(2*q)+ phi_b;
Vfb=phi_m-phi_s;

%solver
dx= 0.1e-9;
NRmax=100;
tol=1e-6;

%mesh 
xox=-tox:dx:0;
xsi= dx:dx:Lsub;
x=[xox,xsi];
N= length (x);
iox= x<=0;
isi=x>0;
i0=find (x==0);

%Newton  Raphson Matrices
ki= k_ox*iox+ k_si*isi;
kip= 0.5*(ki + circshift(ki,-1));
kim= 0.5*(ki + circshift(ki,1));
kipm=kip+kim;
A= diag(kipm,0)- diag (kip(1:end-1),1)- diag (kim(2:end),-1);

 A(1,:)=0; A(1,1)=1;
 b=zeros(N,1);
 
 psi= zeros(N,1);
 
 %fixed oxide charges
 
 rho_ox= zeros(size(xox));
 %rho_ox(round (length(xox)/2))=0*1e12*1e4/dx*q;
 rho_ox(end)= 0*1e12*1e4/dx*q;
 rho_oxfx= rho_ox;
 
 
 %Interface Charges
 %dite=1e14*1e4/q;
 %ditm= 1e12*1e4/q;
 %E0=Eg/(2*log(2*dite/ditm));
 %ditf=@(E) 1.*(dite.*exp(-E/E0)+ dite*exp(-(Eg-E)/E0));
 %input
 ditf=(0*1e12*1e4/Eg);
 
 
 Vgs= (-2.5:50e-3:2.5);
 Vwiggle=10e-3;
 Vright=0;
 sweep_type=1;% 0= LFCV,1=HFCV,2=Fast Sweep
 dit_dpsi=0.1e-3;
 
 Caps=[];
 psi_ss=[];
 
 %Voltage Sweep
 
     for Vg= Vgs
      
     
     
     V1= Vg-Vfb;
     V2=V1+Vwiggle;
     
     Qsub=[];
     
     Vleft=V1;
     for i=1:NRmax
          %figure(1); plot(x,psi);hold on;
          %Generate b
          
          psi_si= psi(isi);
          psi_s= psi(i0);
          E1= q*psi_s+Eg/2-q*phi_b;
          E2=Eg/2;
          rho_dit= (q/dx)*ditf*(E2-E1);
          rho_dit_dpsi=(q/dx)*ditf*(E2-E1-(q*dit_dpsi));
          delrho_dit = 1/dit_dpsi* (rho_dit_dpsi-rho_dit);
         
                   
          p= psub*exp(-psi_si/vt);
          n= nsub*exp(psi_si/vt);
          delp= -1/vt*p;
          deln=1/vt*n;
          
          if sweep_type==2
              minc=0;
              delminc=0;
              ditq=0;
              ditdq=0;
              if sign(Nsub)==1
                  majc=n;
                  delmajc=deln;
              elseif sign (Nsub)==-1
                  majc=p;
                  delmajc=delp;
              end
          else
              ditq=rho_dit;
              ditdq=delrho_dit;
              if sign (Nsub)==1
                  majc=n;
                  delmajc=deln;
                  minc=p;
                  delminc=delp;
              elseif sign(Nsub)==-1
                  majc=p;
                  delmajc=delp;
                  minc=n;
                  delminc=deln;
              end
          end
          rho_ox(end)= rho_oxfx(end)+ditq;
          rho_si= q*(Nsub+sign (Nsub)*(minc-majc));
          
          rho= [rho_ox';rho_si];
          b=rho*dx^2/eps_0;
          b(1)= Vleft;
          b(N)= Vright;
          f= A*psi-b;
          
          %Jacobian Calculation
          delrho_ox= zeros(size(rho_ox));
          delrho_ox(end)= delrho_ox(end)+ ditdq;
          
          delrho_si= q*(sign(Nsub)*(delminc-delmajc));
          
          delrho= [delrho_ox';delrho_si];
          delb=delrho*dx^2/eps_0;
          delb(1)=0;
          delb(N)=0;
          J= A-diag(delb);
          dV= -J\f;
          
          if max(abs(dV))<tol
              break;
          end
          psi=psi+dV;
     end
     
     Q1= sum(rho_si)*dx+ ditq*dx;
  
     minc1=minc;
     ditq1=ditq;
     
     %Wiggle point
     Vleft = V2;
     for i=1:NRmax
         
         %figure(2);plot(x,psi);hold on;
         
         %generate b
         psi_si= psi(isi);
          psi_s= psi(i0);
          E1= q*psi_s+Eg/2-q*phi_b;
          E2=Eg/2;
          rho_dit= (q/dx)*ditf*(E2-E1);
          rho_dit_dpsi=(q/dx)*ditf*(E2-E1-(q*dit_dpsi));
          delrho_dit = 1/dit_dpsi* (rho_dit_dpsi-rho_dit);
          p= psub*exp(-psi_si/vt);
          n= nsub*exp(psi_si/vt);
          delp= -1/vt*p;
          deln=1/vt*n;
          if sweep_type==2
              minc=0;
              delminc=0;
              ditq=0;
              ditdq=0;
              if sign(Nsub)==1
                  majc=n;
                  delmajc=deln;
              elseif sign (Nsub)==-1
                  majc=p;
                  delmajc=delp;
              end
          
          elseif sweep_type==1
              minc=minc1;
              delminc=0;
              ditq=ditq1;
              ditdq=0;
              if sign(Nsub)==1
                  majc=n;
                  delmajc=deln;
              elseif sign (Nsub)==-1
                  majc=p;
                  delmajc=delp;
              end
           elseif sweep_type==0
               ditq=rho_dit;
               ditdq=delrho_dit;
              if sign (Nsub)==1
                  majc=n;
                  delmajc=deln;
                  minc=p;
                  delminc=delp;
              elseif sign(Nsub)==-1
                  majc=p;
                  delmajc=delp;
                  minc=n;
                  delminc=deln;
              end
          end
          rho_ox(end)= rho_oxfx(end)+ditq;
          rho_si= q*(Nsub+sign (Nsub)*(minc-majc));
          
          rho= [rho_ox';rho_si];
          b=rho*dx^2/eps_0;
          b(1)= Vleft;
          b(N)= Vright;
          f= A*psi-b;
          
          %Jacobian Calculation
          delrho_ox= zeros(size(rho_ox));
          delrho_ox(end)= delrho_ox(end) + ditdq;
          delrho_si= q*(sign(Nsub)*(delminc-delmajc));
          
          delrho= [delrho_ox';delrho_si];
          delb=delrho*dx^2/eps_0;
          delb(1)=0;
          delb(N)=0;
          J= A-diag(delb);
          dV= -J\f;
          
          if max(abs(dV))<tol
              break;
          end
          psi=psi+dV;
     end
     
     Q2= sum(rho_si)*dx+ditq*dx;
     
     psi_ss=[psi_ss psi_s];
     Cap= 1/Vwiggle*(Q1-Q2);
     Caps=[Caps Cap];
     end
      
  
 
    figure(4);
    plot(Vgs,Caps*10^2)
    hold on
 

 
 
 
 
 
legend (  'Positive oxide Charge' ,'Zero Oxide Charge');
 xlabel ('Vg (in Volt)');
 ylabel ('C (uF/cm^2)');
 title (' Effect of Fixed Oxide Charge on Interface ');
 

 

 
 
 
 
 
 %legend (  'LFCV' ,'HFCV' ,'Deep Depletion');
 %xlabel ('Vg (in Volt)');
 %ylabel ('C (uF/cm^2)');
 %title ('CV Characteristics');
 
