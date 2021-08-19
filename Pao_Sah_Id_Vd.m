q=1.6e-19;
eps_0=8.85e-12;
kT=26e-3*q;

% semiconductor
Nsub=-3.5e17*1e6;
k_si=12;
ni=1.5e10*1e6;      
Eg=1.1*q;
eps_si=k_si*eps_0;
chi_si=4.05*q;
Na=abs(Nsub);
W=1e-6;
L=1e-6;
muf=200e-4;

% oxide
tox=10e-9;           % oxide thickness
k_ox=4; 
eps_ox=k_ox*eps_0; 
Cox=eps_ox/tox;

phi_m=chi_si/q;
phi_b=-sign(Nsub)*kT/q*log(abs(Nsub)/ni);
phi_s= chi_si/q + Eg/(2*q) + phi_b;
Vfb= phi_m - phi_s;     

% Pao-Sah solution at a Id = f_ps(Vg,Vd)


% Vg=0:0.1:2;
% Vd=[0.1 1];



Vd=0:0.1:5;
Vg=[1 2 3 4 5];

for j=1:length(Vd)
for i=1:length(Vg)
    
psi_s_min= -Vg(i)-Vd(j)-abs(Vfb);
psi_s_max= Vg(i)+Vd(j)+abs(Vfb);
dpsi_s=1e-3;

psi_svec= psi_s_min:dpsi_s:psi_s_max;

IdVs = [];
dV = 1e-2;
del_psi= 10e-3;

for V=0:dV:Vd(j)
    
    f1 = @(psi) ni^2/Na*exp(q*(psi-V)/kT);
    f2 = @(psi) (2*kT*Na/eps_si)^0.5*(q*psi/kT + f1(psi)/Na).^0.5;
    f3 = @(psi_s) Vfb + psi_s + eps_si/Cox*f2(psi_s);
    f1byf2= @(psi) f1(psi)./f2(psi);
    
    Vgs = f3(psi_svec);
    psi_s= interp1(real(Vgs),real(psi_svec),Vg(i));
    
    inn_int=integral(f1byf2, del_psi, psi_s);
    IdVs= [IdVs q*muf*W/L*inn_int];
end

Id(i,j) = sum(IdVs)*dV;
end
end

 
 
%Vg,Vd format
 figure(1);
 
  
 
plot(Vd,Id(5,:))
 hold on 
 plot(Vd,Id(4,:))
 hold on 
 plot(Vd,Id(3,:))
 hold on
 plot(Vd,Id(2,:))
 hold on
 plot(Vd,Id(1,:))
 
 xlabel('Vd');
 ylabel('Id');
title ('Id vs Vd for Pao Sah Model')

legend('Vg=4','Vg=3','Vg=2','Vg=2','Vg=1')
 
 
 