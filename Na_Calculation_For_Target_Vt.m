clear
close all;

% Metal and device parameters
% constants
q=1.6e-19;
eps_0=8.85e-12;
kT=26e-3*q;
Vt=kT/q;


% semiconductor
% Nsub=-4e17*1e6;
k_si=12;
ni=1.5e10*1e6;      
Eg=1.1*q;
eps_si=k_si*eps_0;
chi_si=4.05*q;
% Na=abs(Nsub);
W=1e-6;
L=1e-6;
muf=200e-4;
Vth=0.8;



% oxide
tox=10e-9;           % oxide thickness
k_ox=4; 
eps_ox=k_ox*eps_0; 
Cox=eps_ox/tox;

% Na=4e17*1e6;


Na_init = 4e17*1e6;
Na = Na_init;
tol=1e-6;    
for i = 1:10
    
    F=-Eg/2/q-Vt*log (Na/ni)+2*Vt*log(Na/ni) +((4*eps_si*q* Na* Vt* log (Na/ni))^0.5)/Cox-Vth;
    Fdash=Vt/Na+((4*eps_si*q*Vt)^0.5/Cox)*(0.5*(Na*log(Na/ni))^-0.5*(log(Na/ni)+1));
    del_Na = -F/Fdash;
    if max(abs(del_Na))< tol
    break;
    end
    Na = Na + del_Na ;
end