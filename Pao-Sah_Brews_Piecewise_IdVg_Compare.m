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


Vg=0:0.1:5;
Vd=[0.1 1 2];



% Vd=0:0.1:2;
% Vg=[0.5 1 1.5 2];

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
 
 
 

  
  
  yyaxis left
  semilogy(Vg,Id(:,2))
hold on
yyaxis right
  plot(Vg,Id(:,2))
  hold on
 
  
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
yyaxis left
  semilogy(Vg,Id(:,2))
hold on
yyaxis right
  plot(Vg,Id(:,2))
  hold on
 
 q=1.6e-19;
eps_0=8.85e-12;
kT=26e-3*q;
Vt=26e-3;

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


Vg=0:0.1:5;
Vd=[0.1 1 2];


w=sqrt(2*eps_si*phi_b/(q*Na));
Cd=eps_si/w;
m=1+Cd/Cox;

% Vd=0:0.1:5;
% Vg=1:1:5;

%Threshold voltage
psi_s=2*Vt*log(Na/ni); %psi_s=2*phi_b
psi_ox=(sqrt(2*eps_si*q*Na*psi_s))/Cox; %Voltage drop in oxide
Vth=psi_s+psi_ox+Vfb; %Threshold voltage



for j=1:length(Vd)
for i=1:length(Vg)
    if Vg(i)>=Vth 
    if Vd(j)<=(Vg(i)-Vth)/m
        Id(i,j)=muf*Cox*(W/L)*((Vg(i)-Vth-(m*Vd(j))/2).*Vd(j));
    else
        Id(i,j)= muf*Cox*(W/L)*(Vg(i)-Vth)^2 / (2*m);
        
    end
    else
        Id(i,j)=muf*Cox*(W/L)*(m-1)*(Vt^2)*exp((Vg(i)-Vth)/(m*Vt))*(1-exp(-Vd(j)/Vt));
    end
end
end

 
 

 
yyaxis left
  semilogy(Vg,Id(:,2))
hold on
yyaxis right
  plot(Vg,Id(:,2))
  hold on
 xlabel('Vg');
% ylabel('Id');
%legend ('Pao-Sah','Brews','Piece-Wise');
title ('Comparision:Pao-Sah,Brews and Piece-Wise');
 
 




 
 