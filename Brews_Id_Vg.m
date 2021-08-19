% constants
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

% Brews solution at a Id = f_B(Vg,Vd)


Vg=0:0.1:5;
Vd=[0.1 1 2];



% Vd=0:0.1:2;
% Vg=[0.5 1 1.5 2];

for j=1:length(Vd)
for i=1:length(Vg)
    
QT = @(psi_s) Cox*(Vg(i) - Vfb - psi_s);
QD = @(psi_s) (2*eps_si*q*Na*psi_s).^0.5;
QI = @(psi_s) QT(psi_s) - QD(psi_s);
dVdpsi_s = @(psi_s) 1+2*kT/q*(Cox*QT(psi_s) + eps_si*q*Na)./((QT(psi_s)).^2-(QD(psi_s)).^2);
Vgf = @(psi_s,V) Vfb + psi_s + 1/Cox*(2*eps_si*kT*Na)^0.5*(q*psi_s/kT + ni^2/Na^2*exp(q*(psi_s-V)/kT)).^0.5;


psi_s_min= -Vg(i)-Vd(j)-abs(Vfb);
psi_s_max= Vg(i)+Vd(j)+abs(Vfb);
dpsi_s=1e-3;
psi_svec= psi_s_min:dpsi_s:psi_s_max;

psi_ss = interp1(real(Vgf(psi_svec,0)),real(psi_svec),Vg(i));
psi_sd = interp1(real(Vgf(psi_svec,Vd(j))),real(psi_svec),Vg(i));

intf = @(psi_s) QI(psi_s).*dVdpsi_s(psi_s);

Id(i,j) = muf*W/L*integral(intf,psi_ss,psi_sd);

end
end

figure(1);
 yyaxis left
  semilogy(Vg,Id(:,1))
 hold on 
  semilogy(Vg,Id(:,2))
  hold on
  semilogy(Vg,Id(:,3))
  hold on
  ylabel('Id');
  
  
   yyaxis right
 plot(Vg,Id(:,1))
 hold on 
  plot(Vg,Id(:,2))
 hold on
 plot(Vg,Id(:,3))
 hold on
 
 xlabel('Vg');
 ylabel('Id');
 

 
 title ('Id vs Vg for Brews Model')
 legend('Vd=0.1','Vd=1','Vd=2')
 

