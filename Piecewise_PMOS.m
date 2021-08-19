q=1.6e-19;
eps_0=8.85e-12;
kT=26e-3*q;
Vt=26e-3;

% semiconductor
Nsub=3.5e17*1e6;
k_si=12;
ni=1.5e10*1e6;      
Eg=1.1*q;
eps_si=k_si*eps_0;
chi_si=4.05*q;
Na=abs(Nsub);
W=1e-6;
L=1e-6;
muf=100e-4;

% oxide
tox=10e-9;           % oxide thickness
k_ox=4; 
eps_ox=k_ox*eps_0; 
Cox=eps_ox/tox;

phi_m=chi_si/q+Eg/q;
phi_b=-sign(Nsub)*kT/q*log(abs(Nsub)/ni);
phi_s= chi_si/q + Eg/(2*q) + phi_b;
Vfb= phi_m - phi_s;     

% Pao-Sah solution at a Id = f_ps(Vg,Vd)


Vg=-5:0.1:0;
Vd=[-1 -2];


w=sqrt(2*eps_si*2*abs(phi_b)/(q*Na));
Cd=eps_si/w;
m=1+Cd/Cox;

% Vd=-5:0.001:0;
% Vg=-5:1:-1;


Vth=-0.8;
Vth=abs(Vth);


for j=1:length(Vd)
for i=1:length(Vg)
    if abs(Vg(i))>Vth 
    if abs(Vd(j))<(abs(Vg(i))-Vth)/m
        Id(i,j)=muf*Cox*(W/L)*((abs(Vg(i))-Vth-(m*abs(Vd(j)))/2).*abs(Vd(j)));
    else
        Id(i,j)= muf*Cox*(W/L)*(abs(Vg(i))-Vth)^2 / (2*m);
        
    end
    else
        Id(i,j)=muf*Cox*(W/L)*(m-1)*(Vt^2)*exp((abs(Vg(i))-Vth)/(m*Vt))*(1-exp(-abs(Vd(j))/Vt));
    end
end
end

 
 
%Vg,Vd format
 figure(1);
 
  



plot(Vg,Id(:,1))
hold on
plot(Vg,Id(:,2))
hold on

% plot(Vd,Id(5,:))
%  hold on 
%  plot(Vd,Id(4,:))
%  hold on 
%  plot(Vd,Id(3,:))
%  hold on
%  plot(Vd,Id(2,:))
%  hold on
%  plot(Vd,Id(1,:))
 
 xlabel('Vg');
 ylabel('Id');
title ('Id vs Vg for Piecewise Model-Approximate(PMOS)')
%legend('Vg=-1','Vg=-2','Vg=-3','Vg=-4','Vg=-5')
legend('Vd=-1','Vd=-2')



 
 
 
 
 
 