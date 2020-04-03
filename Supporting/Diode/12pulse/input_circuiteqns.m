function [Acomm,Ncomm,Acond,Ncond,Omegat,Ahat,H] = ...
    input_circuiteqns(h,R,L,Rdc,Ldc,lenHarm,CTFsm,invCTFsm,Tac)

Ncond_coeff = zeros(3,3,12);
Acond = zeros(3,3,12);

Acond(3,3,1:end) = ones(1,1,12)*-(4*R+Rdc)/(4*L+Ldc);
Ncond_coeff(3,3,1:end) = ones(1,1,12)*-1/(4*L+Ldc);

Acond_coeff_1 = [1+2/sqrt(3);-1/sqrt(3);-1-1/sqrt(3)];
Acond_coeff_2 = [1+1/sqrt(3);1/sqrt(3);-1-2/sqrt(3)];
Acond_coeff_3 = [1/sqrt(3);1+1/sqrt(3);-1-2/sqrt(3)];
Acond_coeff_4 = [-1/sqrt(3);1+2/sqrt(3);-1-1/sqrt(3)];
Acond_coeff_5 = [-1-1/sqrt(3);1+2/sqrt(3);-1/sqrt(3)];
Acond_coeff_6 = [-1-2/sqrt(3);1+1/sqrt(3);1/sqrt(3)];
Acond_coeff_7 = [-1-2/sqrt(3);1/sqrt(3);1+1/sqrt(3)];
Acond_coeff_8 = [-1-1/sqrt(3);-1/sqrt(3);1+2/sqrt(3)];
Acond_coeff_9 = [-1/sqrt(3);-1-1/sqrt(3);1+2/sqrt(3)];
Acond_coeff_10 = [1/sqrt(3);-1-2/sqrt(3);1+1/sqrt(3)];
Acond_coeff_11 = [1+1/sqrt(3);-1-2/sqrt(3);1/sqrt(3)];
Acond_coeff_12 = [1+2/sqrt(3);-1-1/sqrt(3);-1/sqrt(3)];

Acond_coeff = [Acond_coeff_1,Acond_coeff_2,Acond_coeff_3,Acond_coeff_4,...
    Acond_coeff_5,Acond_coeff_6,Acond_coeff_7,Acond_coeff_8,...
    Acond_coeff_9,Acond_coeff_10,Acond_coeff_11,Acond_coeff_12];

Ncond_coeff_2_1 = [1+sqrt(3),0,-1];
Ncond_coeff_2_2 = [1,0,-1-sqrt(3)];
Ncond_coeff_2_3 = [0,1,-1-sqrt(3)];
Ncond_coeff_2_4 = [0,1+sqrt(3),-1];
Ncond_coeff_2_5 = [-1,1+sqrt(3),0];
Ncond_coeff_2_6 = [-1-sqrt(3),1,0];
Ncond_coeff_2_7 = [-1-sqrt(3),0,1];
Ncond_coeff_2_8 = [-1,0,1+sqrt(3)];
Ncond_coeff_2_9 = [0,-1,1+sqrt(3)];
Ncond_coeff_2_10 = [0,-1-sqrt(3),1];
Ncond_coeff_2_11 = [1,-1-sqrt(3),0];
Ncond_coeff_2_12 = [1+sqrt(3),-1,0];

Ncond_coeff_2 = [Ncond_coeff_2_1;Ncond_coeff_2_2;Ncond_coeff_2_3;...
    Ncond_coeff_2_4;Ncond_coeff_2_5;Ncond_coeff_2_6;Ncond_coeff_2_7;...
    Ncond_coeff_2_8;Ncond_coeff_2_9;Ncond_coeff_2_10;Ncond_coeff_2_11;...
    Ncond_coeff_2_12];

for i = 1:12
    Acond(1:2,3,i) = CTFsm*Acond_coeff(1:end,i)*Acond(3,3,i);
    Ncond_coeff(1:2,3,i) = CTFsm*Acond_coeff(1:end,i)*Ncond_coeff(3,3,i);
    Ncond_coeff(3,1:2,i) = 1/(4*L+Ldc)*Ncond_coeff_2(i,1:end)*invCTFsm;
    Ncond_coeff(1:2,1:2,i) = CTFsm*Acond_coeff(1:end,i)*Ncond_coeff(3,1:2,i);
end

Acomm = zeros(3,3,12);
Ncomm_coeff = zeros(3,3,12);

Acomm(3,3,1:end)= ones(1,1,12)*-(7/2*R+Rdc)/(7/2*L+Ldc);
Ncomm_coeff(3,3,1:end) = ones(1,1,12)*-1/(7/2*L+Ldc);

Acomm_coeff_1_1 = [1+2/sqrt(3); -1/2-1/sqrt(3); -1/2-1/sqrt(3)];
Acomm_coeff_1_2 = [1+sqrt(3)/2; 0; -1-sqrt(3)/2];
Acomm_coeff_1_3 = [1/2+1/sqrt(3); 1/2+1/sqrt(3); -1-2/sqrt(3)];
Acomm_coeff_1_4 = [0; 1+sqrt(3)/2; -1-sqrt(3)/2];
Acomm_coeff_1_5 = [-1/2-1/sqrt(3); 1+2/sqrt(3); -1/2-1/sqrt(3)];
Acomm_coeff_1_6 = [-1-sqrt(3)/2; 1+sqrt(3)/2; 0];
Acomm_coeff_1_7 = [-1-2/sqrt(3); 1/2+1/sqrt(3); 1/2+1/sqrt(3)];
Acomm_coeff_1_8 = [-1-sqrt(3)/2; 0; 1+sqrt(3)/2];
Acomm_coeff_1_9 = [-1/2-1/sqrt(3); -1/2-1/sqrt(3); 1+2/sqrt(3)];
Acomm_coeff_1_10 = [0; -1-sqrt(3)/2; 1+sqrt(3)/2];
Acomm_coeff_1_11 = [1/2+1/sqrt(3); -1-2/sqrt(3); 1/2+1/sqrt(3)];
Acomm_coeff_1_12 = [1+sqrt(3)/2; -1-sqrt(3)/2; 0];
Acomm_coeff = [Acomm_coeff_1_1,Acomm_coeff_1_2,Acomm_coeff_1_3,...
    Acomm_coeff_1_4,Acomm_coeff_1_5,Acomm_coeff_1_6,Acomm_coeff_1_7,...
    Acomm_coeff_1_8,Acomm_coeff_1_9,Acomm_coeff_1_10,Acomm_coeff_1_11,...
    Acomm_coeff_1_12];

Acomm_coeff_2_1 = [ 0; 1/2+1/sqrt(3); 1/2+1/sqrt(3)];
Acomm_coeff_2_2 = [ -1/3-1/2/sqrt(3); 0; 1/3+1/2/sqrt(3)];
Acomm_coeff_2_3 = -[ 1/2+1/sqrt(3); 1/2+1/sqrt(3); 0];
Acomm_coeff_2_4 = [ 0; -1/3-1/2/sqrt(3); 1/3+1/2/sqrt(3)];
Acomm_coeff_2_5 = [ 1/2+1/sqrt(3); 0; 1/2+1/sqrt(3)];
Acomm_coeff_2_6 = [ 1/3+1/2/sqrt(3); -1/3-1/2/sqrt(3); 0];
Acomm_coeff_2_7 = -[ 0; 1/2+1/sqrt(3); 1/2+1/sqrt(3)];
Acomm_coeff_2_8 = [ 1/3+1/2/sqrt(3); 0; -1/3-1/2/sqrt(3)];
Acomm_coeff_2_9 = [ 1/2+1/sqrt(3); 1/2+1/sqrt(3); 0];
Acomm_coeff_2_10 = [ 0; 1/3+1/2/sqrt(3); -1/3-1/2/sqrt(3)];
Acomm_coeff_2_11 = -[ 1/2+1/sqrt(3); 0; 1/2+1/sqrt(3)];
Acomm_coeff_2_12 = [ -1/3-1/2/sqrt(3); 1/3+1/2/sqrt(3); 0];
Acomm_coeff_2 = -R/L*[Acomm_coeff_2_1,Acomm_coeff_2_2,Acomm_coeff_2_3,...
    Acomm_coeff_2_4,Acomm_coeff_2_5,Acomm_coeff_2_6,Acomm_coeff_2_7,...
    Acomm_coeff_2_8,Acomm_coeff_2_9,Acomm_coeff_2_10,Acomm_coeff_2_11,...
    Acomm_coeff_2_12];

Acomm_coeff_3 = zeros(3,3,12);
Acomm_coeff_3(1:end,1:end,1) = R/L*[ 0,0,0; 0,-1,0; 0,0,-1];
Acomm_coeff_3(1:end,1:end,2) = R/L*[ -1/3 1/3 0; 1/3 -2/3 1/3; 0 1/3 -1/3];
Acomm_coeff_3(1:end,1:end,3) = R/L*[ -1 0 0; 0 -1 0; 0 0 0];
Acomm_coeff_3(1:end,1:end,4) = R/L*[ -2/3 1/3 1/3; 1/3 -1/3 0; 1/3 0 -1/3];
Acomm_coeff_3(1:end,1:end,5) = R/L*[ -1 0 0; 0 0 0; 0 0 -1];
Acomm_coeff_3(1:end,1:end,6) = R/L*[ -1/3 0 1/3; 0 -1/3 1/3; 1/3 1/3 -2/3];
Acomm_coeff_3(1:end,1:end,7) = R/L*[ 0 0 0; 0 -1 0; 0 0 -1];
Acomm_coeff_3(1:end,1:end,8) = R/L*[ -1/3 1/3 0; 1/3 -2/3 1/3; 0 1/3 -1/3];
Acomm_coeff_3(1:end,1:end,9) = R/L*[ -1 0 0; 0 -1 0; 0 0 0];
Acomm_coeff_3(1:end,1:end,10) = R/L*[ -2/3 1/3 1/3; 1/3 -1/3 0; 1/3 0 -1/3];
Acomm_coeff_3(1:end,1:end,11) = R/L*[ -1 0 0; 0 0 0; 0 0 -1];
Acomm_coeff_3(1:end,1:end,12) = R/L*[ -1/3 0 1/3; 0 -1/3 1/3; 1/3 1/3 -2/3];

Ncomm_coeff_2_1 = [ 1+sqrt(3), -1/2, -1/2];
Ncomm_coeff_2_2 = [ 1+sqrt(3), sqrt(3)/2, -1];
Ncomm_coeff_2_3 = [ 1/2, 1/2, -sqrt(3)-1];
Ncomm_coeff_2_4 = [ sqrt(3)/2, 1+sqrt(3), -1];
Ncomm_coeff_2_5 = [ -1/2, 1+sqrt(3), -1/2];
Ncomm_coeff_2_6 = [ -1, 1+sqrt(3), sqrt(3)/2];
Ncomm_coeff_2_7 = [ -1-sqrt(3), 1/2, 1/2];
Ncomm_coeff_2_8 = [ -1, sqrt(3)/2, 1+sqrt(3)];
Ncomm_coeff_2_9 = [ -1/2, -1/2, sqrt(3)+1];
Ncomm_coeff_2_10 = [ sqrt(3)/2, -1, 1+sqrt(3)];
Ncomm_coeff_2_11 = [ 1/2, -1-sqrt(3), 1/2];
Ncomm_coeff_2_12 = [ 1+sqrt(3), -1, sqrt(3)/2];
Ncomm_coeff_2 = [Ncomm_coeff_2_1;Ncomm_coeff_2_2;Ncomm_coeff_2_3;...
    Ncomm_coeff_2_4;Ncomm_coeff_2_5;Ncomm_coeff_2_6;Ncomm_coeff_2_7;...
    Ncomm_coeff_2_8;Ncomm_coeff_2_9;Ncomm_coeff_2_10;Ncomm_coeff_2_11;...
    Ncomm_coeff_2_12];

Ncomm_coeff_3 = zeros(3,3,12);
Ncomm_coeff_3(1:end,1:end,1) = [ 0, 0, 0; 0, 1/2, -1/2; 0, -1/2, 1/2];
Ncomm_coeff_3(1:end,1:end,2) = [ 0, -1/2, 0; 0, 1, 0; 0, -1/2, 0];
Ncomm_coeff_3(1:end,1:end,3) = [ 1/2, -1/2, 0; -1/2, 1/2, 0; 0, 0, 0];
Ncomm_coeff_3(1:end,1:end,4) = [ 1, 0, 0; -1/2, 0, 0; -1/2, 0, 0];
Ncomm_coeff_3(1:end,1:end,5) = [ 1/2, 0, -1/2; 0, 0, 0; -1/2, 0, 1/2];
Ncomm_coeff_3(1:end,1:end,6) = [ 0, 0, -1/2; 0, 0, -1/2; 0, 0, 1];
Ncomm_coeff_3(1:end,1:end,7) = [ 0, 0, 0; 0, 1/2, -1/2; 0, -1/2, 1/2];
Ncomm_coeff_3(1:end,1:end,8) = [ 0, -1/2, 0; 0, 1, 0; 0, -1/2, 0];
Ncomm_coeff_3(1:end,1:end,9) = [ 1/2, -1/2, 0; -1/2, 1/2, 0; 0, 0, 0];
Ncomm_coeff_3(1:end,1:end,10) = [ 1, 0, 0; -1/2, 0, 0; -1/2, 0, 0];
Ncomm_coeff_3(1:end,1:end,11) = [ 1/2, 0, -1/2; 0, 0, 0; -1/2, 0, 1/2];
Ncomm_coeff_3(1:end,1:end,12) = [ 0, 0, -1/2; 0, 0, -1/2; 0, 0, 1];

for i = 1:12
   Ncomm_coeff(3,1:2,i) = 1/(7/2*L+Ldc)*Ncomm_coeff_2(i,1:end)*invCTFsm; 
   Acomm(1:2,3,i) = CTFsm*Acomm_coeff(1:end,i)*Acomm(3,3,i)+...
       CTFsm*Acomm_coeff_2(1:end,i);
   Acomm(1:2,1:2,i) = CTFsm*Acomm_coeff_3(1:end,1:end,i)*invCTFsm;
   Ncomm_coeff(1:2,1:2,i) = 1/L*CTFsm*Ncomm_coeff_3(1:end,1:end,i)*...
       invCTFsm+CTFsm*Acomm_coeff(1:end,i)*Ncomm_coeff(3,1:2,i);
   Ncomm_coeff(1:2,3,i) = CTFsm*Acomm_coeff(1:end,i)*Ncomm_coeff(3,3,i);
end

%**************************************************************************
Ncond = zeros(3,2*lenHarm+2*h+2,12);
for j = 1:12
    for i = 1:2:2*lenHarm
        Ncond(1:end,i:i+1,j) = Ncond_coeff(1:end,1:2,j);
    end
end
%****************************************
for j = 1:12
    for i = 2*lenHarm+1:2:2*lenHarm+lenHarm+1
        Ncond(1:end,i:i+1,j) = [Ncond_coeff(1:end,3,j),zeros(3,1)];
    end
end
%****************************************

Ncomm = zeros(3,3*lenHarm+1,12);

%****************************************
for j = 1:12
    for i = 1:2:2*lenHarm
        Ncomm(1:end,i:i+1,j) = Ncomm_coeff(1:end,1:2,j);
    end
end
%****************************************
%****************************************
for j = 1:12
    for i = 2*lenHarm+1:2:2*lenHarm+lenHarm+1
        Ncomm(1:end,i:i+1,j) = [Ncomm_coeff(1:end,3,j),zeros(3,1)];
    end
end
%****************************************

Omegat = zeros(3*lenHarm+1,3*lenHarm+1);

count = -h;
for i = 1:2:2*lenHarm
    Omegat(i:i+1,i:i+1) = [0 -count;count 0];
    count = count + 1;
end
count = 0;
for i = 2*lenHarm+1:2:2*lenHarm+lenHarm
    Omegat(i:i+1,i:i+1) = [0 -count;count 0];
    count = count + 1;
end

%This loop is used to solve for the alpha that will yield the required
%power sourced by the load.

Ahat = zeros(3+3*lenHarm+1+3*lenHarm+1);

H = zeros(3*lenHarm+1,3);
Mm  = Omegat;

count = 1;
for i = h:-1:-h
    H(count:count+1,1:2) = 1/2/pi*[cos(2*i*pi) sin(2*i*pi);...
        -sin(2*i*pi) cos(2*i*pi)]; %this COULD be wrong, not sure though.
    count = count + 2;
end

for i = count:2:count+2*h+2-1
	H(i,3) = 1/Tac;
end

Ahat(3+3*lenHarm+1+1:end,1:3) = H;
Ahat(4:3+3*lenHarm+1,4:3+3*lenHarm+1) = Omegat;
Ahat(3+3*lenHarm+1+1:end,3+3*lenHarm+1+1:end) = Mm;

