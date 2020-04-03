function [w,Vbase,tol,f,Tac,lenHarm,CTF,CTFsm,...
    invCTF,invCTFsm,Mi,phase,Ts,Vab_coeff,Zi] = ...
    inputs_12pulse_diode(Valphabeta,L,R,Ldc,Rdc,h,Vdc)

A = [1 1 1;1 exp(1i*240*pi/180) exp(1i*120*pi/180);...
    1 exp(1i*120*pi/180) exp(1i*240*pi/180)];
invA = inv(A);

w = 1;     %line frequency in rad/s

Vbase = 1;
tol = 0.001;
f = w/2/pi;
Tac = 1/f;
Vout = Vdc;
V_abc = A*[0;1.0;0.0]; %V_abc = A*[0.1;1+1*1i;0.05+0.05*1i];
Mag_abc = [abs(V_abc(1));abs(V_abc(2));abs(V_abc(3))];
Phase_abc = [angle(V_abc(1));angle(V_abc(2));angle(V_abc(3))];

lenHarm = 2*h+1;
CTF = 2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2)];
CTFsm = CTF(1:2,1:3);
invCTF = inv(CTF);
invCTFsm = invCTF(1:3,1:2);

Mi = zeros(12,3);
for i = 0:1
    Mi(1+i,1:end) = [(-2*i+1)*[0 sqrt(3) 0]*invCTFsm 1];
    Mi(3+i,1:end) = [(-2*i+1)*[-sqrt(3) 0 0]*invCTFsm 1];
    Mi(5+i,1:end) = [(-2*i+1)*[0 0 sqrt(3)]*invCTFsm 1];
    Mi(7+i,1:end) = [(-2*i+1)*[0 -sqrt(3) 0]*invCTFsm 1];
    Mi(9+i,1:end) = [(-2*i+1)*[sqrt(3) 0 0]*invCTFsm 1];
    Mi(11+i,1:end) = [(-2*i+1)*[0 0 -sqrt(3)]*invCTFsm 1];
end

Vab_coeff1 = [-L/(4*L+Ldc)*(1+sqrt(3)),-1+L/(4*L+Ldc),1]*invCTFsm;
Vab_coeff2 = [-L/(4*L+Ldc)*(1+sqrt(3)),-sqrt(3),L/(4*L+Ldc)]*invCTFsm;
Vab_coeff3 = [1-L/(4*L+Ldc),-1,L/(4*L+Ldc)*(1+sqrt(3))]*invCTFsm;
Vab_coeff4 = [sqrt(3)-L/(4*L+Ldc)*sqrt(3),-L/(4*L+Ldc)*(1+sqrt(3)),...
    L/(4*L+Ldc)]*invCTFsm;
Vab_coeff5 = [1,-L/(4*L+Ldc)*(1+sqrt(3)),L/(4*L+Ldc)-1]*invCTFsm;
Vab_coeff6 = [L/(4*L+Ldc),-(1+sqrt(3))*L/(4*L+Ldc),-sqrt(3)]*invCTFsm;
Vab_coeff7 = [L/(4*L+Ldc)*(1+sqrt(3)),-L/(4*L+Ldc)+1,-1]*invCTFsm;
Vab_coeff8 = [L/(4*L+Ldc),sqrt(3)-L/(4*L+Ldc)*sqrt(3),...
    -L/(4*L+Ldc)*(1+sqrt(3))]*invCTFsm;
Vab_coeff9 = [L/(4*L+Ldc)-1,1,-L/(4*L+Ldc)*(1+sqrt(3))]*invCTFsm;
Vab_coeff10 = [-sqrt(3),L/(4*L+Ldc),-L/(4*L+Ldc)*(1+sqrt(3))]*invCTFsm;
Vab_coeff11 = [-1,L/(4*L+Ldc)*(1+sqrt(3)),-L/(4*L+Ldc)+1]*invCTFsm;
Vab_coeff12 = [-L/(4*L+Ldc)*(1+sqrt(3)),L/(4*L+Ldc),...
    -L/(4*L+Ldc)*sqrt(3)+sqrt(3)]*invCTFsm;

Vab_coeff = [Vab_coeff1;Vab_coeff2;Vab_coeff3;Vab_coeff4;Vab_coeff5;...
    Vab_coeff6;Vab_coeff7;Vab_coeff8;Vab_coeff9;Vab_coeff10;Vab_coeff11;...
    Vab_coeff12];

%end initializations************************

Vsabc_phasor = zeros(3,1);
for i = 1:3
    [x y] = pol2cart(Phase_abc(i),Mag_abc(i));
    Vsabc_phasor(i) = x + 1i*y;
end

temp = [1 -1 0;0 1 -1;-1 0 1]*[Vsabc_phasor(1);Vsabc_phasor(2);Vsabc_phasor(3)];

%Testing out new algorithm for switching times**********************
symMatrix = invA*Vsabc_phasor;

posSeq = symMatrix(2);
negSeq = symMatrix(3);

phase = angle(Valphabeta(2*h+3)+1i*Valphabeta(2*h+4));
Ts = zeros(24,1);
Ts(1:end) = [1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12];

Zi = zeros(3*lenHarm+1,1);

Zi(1:2*lenHarm,1) = Valphabeta;
Vdc_vect_dc = Vout;%0.5;

%dc comps
Zi(2*lenHarm+1) = Vdc_vect_dc;

end