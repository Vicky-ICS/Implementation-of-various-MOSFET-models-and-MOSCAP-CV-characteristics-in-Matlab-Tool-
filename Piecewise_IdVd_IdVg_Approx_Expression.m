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


% Vg=0:0.1:5;
% Vd=[0.1 1 ];


w=sqrt(2*eps_si*2*phi_b/(q*Na));
Cd=eps_si/w;
m=1+Cd/Cox;

Vd=0:0.1:5;
Vg=1:1:5;

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

 
 
%Vg,Vd format
 figure(1);
 
  
% yyaxis left
% semilogy(Vg,Id(:,1))
% hold on 
% semilogy(Vg,Id(:,2))
% hold on
% 
% yyaxis right 
% plot(Vg,Id(:,1))
% hold on
% plot(Vg,Id(:,2))
% hold on

plot(Vd,Id(5,:))
 hold on 
 plot(Vd,Id(4,:))
 hold on 
 plot(Vd,Id(3,:))
 hold on
 plot(Vd,Id(2,:))
 hold on
 plot(Vd,Id(1,:))
 
 %xlabel('Vd');
 %ylabel('Id');
%title ('Id vs Vg for Piecewise Model-Approximate ')
%legend('Vg=5','Vg=4','Vg=3','Vg=2','Vg=1')
%legend('Vd=0.1','Vd=1','Vd=2')



 
 
 
 
 
 