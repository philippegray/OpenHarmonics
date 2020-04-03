function [Ia Ib Ic fundCurrents phaseAng gamma mu Pload] =...
    solveDiode6(Valphabeta,L,R,Rdc,Ldc,h,Vdc)

global SymMtx invSymMtx

%Initializations***************************
countIter = 0;
w = 1;     %line frequency in rad/s
Rout = Rdc;
Lout = Ldc; %at 100 or even 10 for that matter, the current just keeps rising unbounded it seems.
tolerance = 0.001;
f = w/2/pi;
Tac = 1/f;
Vout = Vdc; %assuming first that Vout is actually just a DC value...it won't be necessarily though so I will change this later.  This problem is directly analogous with the Idc FCM issues encountered in the VSC model.
lenHarm = 2*h+1;
Ra = R;
Rb = R;
Rc = R;
La = L;
Lb = L;
Lc = L;

V_abc = SymMtx*[0;Valphabeta(lenHarm+2)+1i*Valphabeta(lenHarm+3);...
    Valphabeta(2*h-1)+1i*Valphabeta(2*h)];
Mag_abc = [abs(V_abc(1));abs(V_abc(2));abs(V_abc(3))];
Phase_abc = [angle(V_abc(1));angle(V_abc(2));angle(V_abc(3))];

CTF = 2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2)];
invCTF = inv(CTF);
invCTFsmall = invCTF(1:3,1:2);
arrPosArr = [1;2;3;1;2;3];

Mi = zeros(6,3);
Mi(2,1:end) = [1 0 0];
Mi(3,1:end) = -1*[-1/2 -sqrt(3)/2 0];
Mi(1,1:end) = -1*[-1/2 sqrt(3)/2 0];
Mi(5,1:end) = -1*[1 0 0];
Mi(6,1:end) = [-1/2 -sqrt(3)/2 0];
Mi(4,1:end) = [-1/2 sqrt(3)/2 0];
Mi_gamma = [0 0 1];
%end initializations************************

Vsabc_phasor = zeros(3,1);
for i = 1:3
    [x y] = pol2cart(Phase_abc(i),Mag_abc(i));
    Vsabc_phasor(i) = x + 1i*y;
end

temp = [1 -1 0;0 1 -1;-1 0 1]*[Vsabc_phasor(1);Vsabc_phasor(2);Vsabc_phasor(3)];

%Testing out new algorithm for switching times**********************
symMatrix = invSymMtx*Vsabc_phasor;

posSeq = symMatrix(2);

phaseAng = angle(posSeq)-0;%-30/180*pi;%150/180*pi;
Ts = zeros(18,1);
Ts(1:end) = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6];

Zi = zeros(3*lenHarm+1,1);

Zi(1:2*lenHarm) = Valphabeta;
Zi(2*lenHarm+1) = Vdc;

%Now that the Ts(1,1) and alpha are known and the input vector has been
%provided...we can use a rotation matrix to move the input vector to a time
%equal to Ts(1,1) + alpha

%Ts(*,*) contains the swithing times and the 2 phase voltage conduction for
%each conduction interval from 0 -> 2*pi
%An example: Ts(i,1) contains the next switching time, and Ts(i,2) contains
%the 2 phase voltage conduction for Ts(i-1,1) --> Ts(i,1)

Acomm = zeros(3,3,6);
Acond = zeros(3,3,6);
Ncond = zeros(3,3*lenHarm+1,6);
%Transition from Vab -> Vac
Acomm(1:2,1:end,1) = [-(3/2*R+Rout)/(3/2*L+Lout),0,0;...
    0,-R/L,0];
Acomm(3,1:end,1) = [0,0,-(3/2*R+Rout)/(3/2*L+Lout)];
%Transition from Vca -> Vcb
Acomm(1:2,1:end,5) = [1/(2*L*L+4*(L+Lout)*L)*(-2*L*R-3*(L+Lout)*R-L*(R+Rout)),...
    1/(2*L*L+4*(L+Lout)*L)*(sqrt(3)*(L+Lout)*R-sqrt(3)*L*(R+Rout)),0;...
    1/(2*sqrt(3)*L*L+4*sqrt(3)*L*(L+Lout))*(3*(L+Lout)*R-3*L*(R+Rout)),...
    1/(2*sqrt(3)*L*L+4*sqrt(3)*L*(L+Lout))*(-2*sqrt(3)*L*R-sqrt(3)*(L+Lout)*R-3*sqrt(3)*L*(R+Rout)),0];
Acomm(3,1:end,5) = [0,0,-(3/2*R+Rout)/(3/2*L+Lout)];
%Transition from Vbc -> Vba
Acomm(1:2,1:end,3) = [-(2*L*R+3*(L+Lout)*R+L*(R+Rout))/(2*L*L+4*(L+Lout)*L),...
    -(sqrt(3)*(L+Lout)*R-sqrt(3)*L*(R+Rout))/(2*L*L+4*(L+Lout)*L),0;...
    -(3*(L+Lout)*R-3*L*(R+Rout))/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L),...
    (-2*sqrt(3)*L*R-sqrt(3)*(L+Lout)*R-sqrt(3)*3*L*(R+Rout))/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L),0];
Acomm(3,1:end,3) = [0,0,-(3/2*R+Rout)/(3/2*L+Lout)];
%Transition from Vba -> Vca
Acomm(1:end,1:end,4) = Acomm(1:end,1:end,1);
%Transition from Vac -> Vbc 
Acomm(1:end,1:end,2) = Acomm(1:end,1:end,5);
%Transition from Vcb -> Vab
Acomm(1:end,1:end,6) = Acomm(1:end,1:end,3);

%Now initializing the conduction A arrays

Acond(1:2,1:end,6) = [-(Ra+Rb+Rout)/(La+Lb+Lout),0,0;...
    0,-(Ra+Rb+Rout)/(La+Lb+Lout),0];
Acond(3,1:end,6) = [0,0,-(Ra+Rb+Rout)/(La+Lb+Lout)];
Acond(1:2,1:end,2) = [-(Ra+Rb+Rout)/(La+Lb+Lout),0,0;...
    0,-(Ra+Rb+Rout)/(La+Lb+Lout),0];
Acond(3,1:end,2) = [0,0,-(Rc+Rb+Rout)/(Lc+Lb+Lout)];
Acond(1:2,1:end,4) = [-(Ra+Rb+Rout)/(La+Lb+Lout),0,0;...
    0,-(Ra+Rb+Rout)/(La+Lb+Lout),0];
Acond(3,1:end,4) = [0,0,-(Ra+Rc+Rout)/(La+Lc+Lout)];
Acond(1:end,1:end,1) = Acond(1:end,1:end,4);
Acond(1:end,1:end,3) = Acond(1:end,1:end,6);
Acond(1:end,1:end,5) = Acond(1:end,1:end,2);
%****************************************
%first row is for the disalpha/dt, column 1 vsalpha, column 2 vsbeta,
%column 3 v_{dc}
%second row is for the disbeta/dt
%third dimension is for differentiating between the different intervals. -
%CHECK THAT THIS IS INDEED RIGHT I COULD HAVE FORGOTTEN AND DONE THIS WITH
%RESPECT TO THE FIRST SECOND THIRD INTERVAL WHEN IT SHOULD HAVE BEEN WITH
%RESPECT TO THE NUMBER OF THE 
coeffArr = zeros(3,3,6);
coeffArr(1,1,6) = 3/2/(La+Lb+Lout);
coeffArr(1,2,6) = -sqrt(3)/2/(La+Lb+Lout);
coeffArr(1,3,6) = -1/(La+Lb+Lout);
coeffArr(2,1,6) = -sqrt(3)/2/(La+Lb+Lout);
coeffArr(2,2,6) = 1/2/(La+Lb+Lout);
coeffArr(2,3,6) = 1/sqrt(3)/(La+Lb+Lout);
coeffArr(3,1,6) = 3/2/(La+Lb+Lout);
coeffArr(3,2,6) = -sqrt(3)/2/(La+Lb+Lout);
coeffArr(3,3,6) = -1/(La+Lb+Lout);
coeffArr(2,2,2) = 2/(Lb+Lc+Lout);
coeffArr(2,3,2) = -2/sqrt(3)/(Lb+Lc+Lout);
coeffArr(3,2,2) = sqrt(3)/(Lc+Lb+Lout);
coeffArr(3,3,2) = -1/(Lc+Lb+Lout);
coeffArr(1,1,4) = 3/2/(La+Lc+Lout);
coeffArr(1,2,4) = sqrt(3)/2/(Lc+La+Lout);
coeffArr(1,3,4) = 1/(La+Lc+Lout);
coeffArr(2,1,4) = sqrt(3)/2/(La+Lc+Lout);
coeffArr(2,2,4) = 1/2/(La+Lc+Lout);
coeffArr(2,3,4) = 1/sqrt(3)/(La+Lc+Lout);
coeffArr(3,1,4) = -3/2/(La+Lc+Lout);
coeffArr(3,2,4) = -sqrt(3)/2/(La+Lc+Lout);
coeffArr(3,3,4) = -1/(La+Lb+Lout);
coeffArr(1:end,1:end,1) = coeffArr(1:end,1:end,4);
coeffArr(1:end,1:end,3) = coeffArr(1:end,1:end,6);
coeffArr(1:end,1:end,5) = coeffArr(1:end,1:end,2);
coeffArr(1:2,3,1) = coeffArr(1:2,3,1)*-1;
coeffArr(1:2,3,3) = coeffArr(1:2,3,3)*-1;
coeffArr(1:2,3,5) = coeffArr(1:2,3,5)*-1;
coeffArr(3,1,1) = coeffArr(3,1,1)*-1;
coeffArr(3,2,1) = coeffArr(3,2,1)*-1;
coeffArr(3,1,3) = coeffArr(3,1,3)*-1;
coeffArr(3,2,3) = coeffArr(3,2,3)*-1;
coeffArr(3,1,5) = coeffArr(3,1,5)*-1;
coeffArr(3,2,5) = coeffArr(3,2,5)*-1;
%****************************************

for j = 1:6
    for i = 1:2:2*lenHarm
        Ncond(1:end,i:i+1,j) = coeffArr(1:end,1:2,j);
    end
end
%****************************************
for j = 1:6
    for i = 2*lenHarm+1:2:2*lenHarm+lenHarm+1
        Ncond(1:end,i:i+1,j) = [coeffArr(1,3,j),0;coeffArr(2,3,j),0;coeffArr(3,3,j),0];
    end
end
%****************************************

coeffArr_comm = zeros(3,3,6);
Ncomm = zeros(3,3*lenHarm+1,6);

coeffArr_comm(1:end,1:end,1) = [3/2/(3/2*L+Lout) 0 -1/(3/2*L+Lout);...
    0 1/L 0;3/2/(3/2*L+Lout) 0 -1/(3/2*L+Lout)];
coeffArr_comm(1:end,1:end,5) = [(3*L+3*(L+Lout))/(2*L*L+4*(L+Lout)*L),...
    -sqrt(3)*Lout/(2*L*L+4*(L+Lout)*L),2*L/(2*L*L+4*(L+Lout)*L);...
    -3*Lout/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L),...
    (5*sqrt(3)*L+sqrt(3)*(L+Lout))/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L),...
    6*L/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L);-3/4/(3/2*L+Lout),...
    -3*sqrt(3)/4/(3/2*L+Lout),-1/(3/2*L+Lout)];
coeffArr_comm(1:end,1:end,3) = [(3*L+3*(L+Lout))/(2*L*L+4*(L+Lout)*L),...
    sqrt(3)*Lout/(2*L*L+4*(L+Lout)*L),2*L/(2*L*L+4*(L+Lout)*L);...
    (3*(L+Lout)-3*L)/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L),...
    (5*sqrt(3)*L+sqrt(3)*(L+Lout))/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L),...
    -6*L/(2*sqrt(3)*L*L+4*sqrt(3)*(L+Lout)*L);-3/4/(3/2*L+Lout),...
    3*sqrt(3)/4/(3/2*L+Lout) -1/(3/2*L+Lout)];
coeffArr_comm(1:end,1:end,2) = coeffArr_comm(1:end,1:end,5);
coeffArr_comm(1:end,1:end,4) = coeffArr_comm(1:end,1:end,1);
coeffArr_comm(1:end,1:end,6) = coeffArr_comm(1:end,1:end,3);
coeffArr_comm(1:2,3,2) = coeffArr_comm(1:2,3,2)*-1;
coeffArr_comm(1:2,3,4) = coeffArr_comm(1:2,3,4)*-1;
coeffArr_comm(1:2,3,6) = coeffArr_comm(1:2,3,6)*-1;
coeffArr_comm(3,1:2,2) = coeffArr_comm(3,1:2,2)*-1;
coeffArr_comm(3,1:2,4) = coeffArr_comm(3,1:2,4)*-1;
coeffArr_comm(3,1:2,6) = coeffArr_comm(3,1:2,6)*-1;

%****************************************
for j = 1:6
    for i = 1:2:2*lenHarm
        Ncomm(1:end,i:i+1,j) = coeffArr_comm(1:3,1:2,j);
    end
end
%****************************************
%****************************************
for j = 1:6
    for i = 2*lenHarm+1:2:2*lenHarm+lenHarm+1
        Ncomm(1:end,i:i+1,j) = [coeffArr_comm(1,3,j),0;...
            coeffArr_comm(2,3,j),0;coeffArr_comm(3,3,j),0];
    end
end
%****************************************

%I'm assuming everything previous to this point is correct...I'm pretty
%sure it is as I have checked.  That being said,there may be an error
%somewhere in one of the coefficient arrays.

%mu = pi/6*ones(6,1); %Initialization
mu = (pi/3-0)*ones(6,1);%(pi/3-0.2)*ones(6,1);
gamma = pi/3*ones(6,1);
phaseAng = 0.501;
%The following for loop calculates the delta mu(1)...mu(3) terms in the
%Jacobian matrix formulation.
At = zeros(3,3);
Nt = zeros(3,3*lenHarm+1); %why is this 4
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

H = zeros(3*lenHarm+1,3);
Mm  = Omegat;

count = 1;
for i = h:-1:-h
    H(count:count+1,1:2) = 1/2/pi*[cos(i*2*pi) sin(i*2*pi);-sin(i*2*pi) cos(i*2*pi)];
    count = count + 2;
end

for i = count:2:count+2*h+2-1
	H(i,3) = 1/Tac;
end

Mcond = zeros(3+3*lenHarm+1,3+3*lenHarm+1);
Mcomm = zeros(3+3*lenHarm+1,3+3*lenHarm+1);
Mt_arr = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6);
dMp_dmu = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6,6); %3rd column contains the partial
dMp_dgamma = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6,6); %3rd column contains the partial
Mp = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6);
Ap = zeros(3,3,6);
Np = zeros(3,3*lenHarm+1,6);
Omegap = zeros(3*lenHarm+1,3*lenHarm+1,6);

Ahat = zeros(3+3*lenHarm+1+3*lenHarm+1);

Ahat(3+3*lenHarm+1+1:end,1:3) = H;
Ahat(4:3+3*lenHarm+1,4:3+3*lenHarm+1) = Omegat;
Ahat(3+3*lenHarm+1+1:end,3+3*lenHarm+1+1:end) = Mm;
phase_init = phaseAng;

% for u = -100:1:900
%     phase = phase_init+u/1000*2*pi;
%phase = phase+0.0944;
while(true)
    
    %This loop is used to solve for mu(1) --> mu(6) & gamma(1) -> gamma(6)
    %through a Newton Raphson Iterative solving technique.
    
    RotationMatrix = zeros(3*lenHarm+1,3*lenHarm+1,7);
    Zi_arr = zeros(3*lenHarm+1,7);
    gamma_sum = 0;
    for k = 1:7
        count = 1;
        for i  = -h:1:h
            RotationMatrix(count:count+1,count:count+1,k) = ...
                [cos(i*(gamma_sum-phaseAng)),-sin(i*(gamma_sum-phaseAng));...
                sin(i*(gamma_sum-phaseAng)),cos(i*(gamma_sum-phaseAng))];
            count = count + 2;
        end
        for i = 0:h
            RotationMatrix(count:count+1,count:count+1,k) = ...
                [cos(i*(gamma_sum-phaseAng)),-sin(i*(gamma_sum-phaseAng));...
                sin(i*(gamma_sum-phaseAng)),cos(i*(gamma_sum-phaseAng))];
            count = count + 2;
        end
        if k < 7
            gamma_sum = gamma_sum + gamma(k);
        end
        Zi_arr(1:end,k) = RotationMatrix(1:end,1:end,k)*Zi;
    end
    
    for j = 1:6
        Mp(1:end,1:end,j) = eye(length(Mp));
        for k = 1:6 %k = 1->6 for mu, k = 7->12 for gamma
            dMp_dmu(1:end,1:end,j,k) = eye(3+3*lenHarm+1);%zeros(length(derivMp));%
            dMp_dgamma(1:end,1:end,j,k) = eye(3+3*lenHarm+1);%zeros(length(derivMp));%
        end
    end
    
    expmCommMu = zeros(length(Ahat),length(Ahat),6);
    expmCondMu = zeros(length(Ahat),length(Ahat),6);
    expmCondMuPi_3 = zeros(length(Ahat),length(Ahat),6);
    
    Phi = eye(3+3*lenHarm+1+3*lenHarm+1);
    for i = 1:6
        Ahat(1:3,1:3+3*lenHarm+1) = [Acomm(1:end,1:end,i),Ncomm(1:end,1:end,i)];
        expmCommMu(1:end,1:end,i) = expm(Ahat*mu(i));
        Ahat(1:3,1:3+3*lenHarm+1) = [Acond(1:end,1:end,i),Ncond(1:end,1:end,i)];
        expmCondMu(1:end,1:end,i) = expm(-Ahat*mu(i));
        expmCondMuPi_3(1:end,1:end,i) = expm(Ahat*gamma(i));
        
        Phi = expmCondMuPi_3(1:end,1:end,i)*...
            expmCondMu(1:end,1:end,i)*expmCommMu(1:end,1:end,i)*Phi;
    end
    
    maxPos = 3+3*lenHarm+1;
    
    for i = 1:6
        for j = 0:5
            Mcond(1:3,1:3) = Acond(1:end,1:end,Ts(i+j));
            Mcond(1:3,4:end) = Ncond(1:end,1:end,Ts(i+j));
            Mcond(4:end,4:end) = Omegat;
            Mcomm(1:3,1:3) = Acomm(1:end,1:end,Ts(i+j));
            Mcomm(1:3,4:end) = Ncomm(1:end,1:end,Ts(i+j));
            Mcomm(4:end,4:end) = Omegat;
            
            Mp(1:end,1:end,i) = expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*Mp(1:end,1:end,i);
            
            for k = 1:6
                if k == (j+1)
                    dMp_dmu(1:end,1:end,k,i) = ...
                        expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                        (-Mcond*expmCondMu(1:maxPos,1:maxPos,Ts(i+j))+...
                        expmCondMu(1:maxPos,1:maxPos,Ts(i+j))*Mcomm)*...
                        expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                        dMp_dmu(1:end,1:end,k,i);
                    
                    dMp_dgamma(1:end,1:end,k,i) = Mcond*...
                        expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                        dMp_dgamma(1:end,1:end,k,i);
                else
                    dMp_dmu(1:end,1:end,k,i) = ...
                        expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                        dMp_dmu(1:end,1:end,k,i);
                    
                    dMp_dgamma(1:end,1:end,k,i) = ...
                        expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                        dMp_dgamma(1:end,1:end,k,i);
                end
            end
        end
    end
    
    %the 4th dimension is whether it be M1,M2,M3 (each has an offset to the initial vector)
    %the 3rd dimension is for the derivative with respect to mu1,mu2,mu3
    dMpgamma = zeros(length(Mp),length(Mp),6,6);
    dMpmu = zeros(length(Mp),length(Mp),6,6);
    
    for i = 0:5
        for j = 1:6
            dMpgamma(1:end,1:end,Ts(j+i),i+1) = dMp_dgamma(1:end,1:end,j,i+1);
            dMpmu(1:end,1:end,Ts(j+i),i+1) = dMp_dmu(1:end,1:end,j,i+1);
        end
    end
    
    %the Mp derivative terms should be fine now.
    for i = 1:6
        Ap(1:end,1:end,i) = Mp(1:3,1:3,i);
        Np(1:end,1:end,i) = Mp(1:3,4:end,i);
        Omegap(1:end,1:end,i) = Mp(4:end,4:end,i);
    end
    
    derivMt = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6);
    
    %M2 SOLVE
    %This block calculates the derivative of the Mt matrix elements.
    
    for i = 1:6
        Mcomm(1:3,1:3) = Acomm(1:end,1:end,i);
        Mcomm(1:3,4:end) = Ncomm(1:end,1:end,i);
        Mcomm(4:end,4:end) = Omegat;
        
        Mt_arr(1:end,1:end,i) = expmCommMu(1:maxPos,1:maxPos,i);
        derivMt(1:end,1:end,i) = Mcomm*Mt_arr(1:end,1:end,i);
    end
    
    %Derivation of FCM
    Hp = Phi(4+3*lenHarm+1:end,1:3);
    Ap_mat = Phi(1:3,1:3);
    Np_mat = Phi(1:3,4:3*lenHarm+4);
    Qp = Phi(3*lenHarm+4+1:end,4:3*lenHarm+4);
    FCM = Hp*(inv(eye(3,3)-Ap_mat))*Np_mat+Qp;
    %testing new FCM formulation since currently, the third column of Ap
    %contains all 0's --> which implies a non-invertable matrix.
    
    Ahat2 = Ahat;
    
    I = FCM*Zi_arr(1:end,1);
    
    %Is Ialpha + Ibeta - the calculated value similar!!!...if not then redo
    %loop
    
    derivRotation = zeros(length(FCM),length(FCM),7);
    gamma_sum = 0;
    for j = 1:7
        count = 1;
        for i  = -h:1:h
            derivRotation(count:count+1,count:count+1,j) = [-i*sin(i*(gamma_sum-phaseAng)),...
                -i*cos(i*(gamma_sum-phaseAng));i*cos(i*(gamma_sum-phaseAng)) -i*sin(i*(gamma_sum-phaseAng))];
            count = count + 2;
        end
        for i  = 0:1:h
            derivRotation(count:count+1,count:count+1,j) = [-i*sin(i*(gamma_sum-phaseAng)),...
                -i*cos(i*(gamma_sum-phaseAng));i*cos(i*(gamma_sum-phaseAng)) -i*sin(i*(gamma_sum-phaseAng))];
            count = count + 2;
        end
        if j < 7
            gamma_sum = gamma_sum + gamma(j);
        end
    end
    
    Idc = I(2*lenHarm+1);
    Vdc = Zi(2*lenHarm+1);
    
    %I also need the correct Mcomm...
    %There has to be a conditional if statement that if j ==
    %arrPosArr(Ts(6,2)), then for M1 M2 M3 the values have to be saved for
    %that derivative term.  Actually, might not even need to use a loop for
    %this...could just directly add the terms as no further derivateves are
    %required.
    %**************************************************
    
    %This needs a little fix up!*******************************************
    J11 = zeros(6,6);
    J12 = zeros(6,7);
    J21 = zeros(7,6);
    J22 = zeros(7,7);
    %there is going to be four regions
    
    %******************
    %At this point the Matrices should be solved for and I should be able
    %to immediatly obtain M1, M2, M3
    
    M = zeros(13,1);
    for i = 1:6
        At = Mt_arr(1:3,1:3,i);
        Nt = Mt_arr(1:3,4:end,i);
        M(i) = 0-Mi(i,1:end)*(At*inv(eye(3,3)-Ap(1:3,1:3,i))*...
            Np(1:3,1:end,i)+Nt)*Zi_arr(1:end,i);
    end
    
    Vabc_coeff = zeros(6,2);
    Vabc_coeff(2,1:end) = [1-L/(2*L+Lout),-1,L/(2*L+Lout)]*invCTFsmall;
    Vabc_coeff(3,1:end) = [1,-L/(2*L+Lout),-1+L/(2*L+Lout)]*invCTFsmall;
    Vabc_coeff(4,1:end) = [L/(2*L+Lout),1-L/(2*L+Lout),-1]*invCTFsmall;
    Vabc_coeff(5,1:end) = [-1+L/(2*L+Lout),1,-L/(2*L+Lout)]*invCTFsmall;
    Vabc_coeff(6,1:end) = [-1,L/(2*L+Lout),1-L/(2*L+Lout)]*invCTFsmall;
    Vabc_coeff(1,1:end) = [-L/(2*L+Lout),-1+L/(2*L+Lout),1]*invCTFsmall;
    %Vabc_coeff(1,1:end) = -1/(2*L+Lout)*[L,-(L+Lout) -2*L]*invCTFsmall;
    constMtx = zeros(3,3*lenHarm+1);
    for i = 1:2
        for j = i:2:2*lenHarm
            constMtx(i,j) = 1;
        end
    end
    for i = 1:2:lenHarm+1
        constMtx(3,2*lenHarm+i) = 1;
    end
    
    for i = 1:7 %need to use Ts b/c I am solving for 7 intervals > num of intervals/period
        Idc_temp = Mi_gamma*(inv(eye(3,3)-Ap(1:3,1:3,Ts(i)))*...
            Np(1:3,1:end,Ts(i)))*Zi_arr(1:end,i);
        
%   ORIGINAL       M(i+6) = 0 - ((L*(2*R+Rout)/(2*L+Lout)-R)*Idc_temp+...
%              [Vabc_coeff(Ts(i),1:end),L/(2*L+Lout)]*constMtx*Zi_arr(1:end,i));
         M(i+6) = 0 - ((L*Rout-R*Lout)/(2*L+Ldc)*Idc_temp+...
             [Vabc_coeff(Ts(i),1:end),L/(2*L+Lout)]*constMtx*Zi_arr(1:end,i));
    end %why does Zi_arr have to be Ts(i+1) - shouldn't it be just i..
    %calculating M(13)*********
    
    %**************************
    %Calculation of J******************************************************
    for i = 1:6
        for j = 1:6 
            %Entries in Quad. (1,1)
            J11(i,j) = calcJ11(Ap(1:end,1:end,i),Np(1:end,1:end,i),...
                Mt_arr(1:end,1:end,i),derivMt(1:end,1:end,i),...
                dMpmu(1:end,1:end,j,i),Zi_arr(1:end,i),Mi(i,1:end),i-j,lenHarm);
            if j < i
                condMtx = derivRotation(1:end,1:end,i);
            else
                condMtx = zeros(3*lenHarm+1,3*lenHarm+1);
            end
            %Entries in Quad. (1,2)
            J12(i,j) = calcJ12(Ap(1:end,1:end,i),Np(1:end,1:end,i),...
                Mt_arr(1:end,1:end,i),dMpgamma(1:end,1:end,j,i),...
                Zi_arr(1:end,i),Zi,Mi(i,1:end),condMtx);
        end
        condMtx = derivRotation(1:end,1:end,i);
        J12(i,7) = calcJ13(Mt_arr(1:end,1:end,i),Ap(1:end,1:end,i),Np(1:end,1:end,i),...
            Zi,Mi(i,1:end),condMtx);
    end
    for i = 1:7
        for j = 1:6 
            %Entries in Quad. (2,1)
            J21(i,j) = calcJ21(Ap(1:end,1:end,Ts(i)),Np(1:end,1:end,Ts(i)),... %might be looking @wrong value??
                dMpmu(1:end,1:end,j,Ts(i)),Zi_arr(1:end,i),Mi_gamma,L,R,Rout,Lout);
            if j < i
                condMtx = derivRotation(1:end,1:end,i);
            else
                condMtx = zeros(3*lenHarm+1,3*lenHarm+1);
            end
            J22(i,j) = calcJ22(Ap(1:end,1:end,Ts(i)),Np(1:end,1:end,Ts(i)),...
                dMpgamma(1:end,1:end,j,Ts(i)),Zi_arr(1:end,i),Zi,Mi_gamma,...
                condMtx,[Vabc_coeff(Ts(i),1:end),L/(2*L+Lout)]*constMtx,Lout,Rout,L,R);
        end
        condMtx = derivRotation(1:end,1:end,i);
        J22(i,7) = calcJ23(Ap(1:end,1:end,Ts(i)),Np(1:end,1:end,Ts(i)),...
            Zi,Mi_gamma,condMtx,[Vabc_coeff(Ts(i),1:end),L/(2*L+Lout)]*constMtx,Lout,Rout,L,R);
    end
    
    J = [J11,J12;J21,J22];
    %**********************************************************************
    
    %mu = mu + J1\M(1:6);
    %gamma = gamma + J2\M(7:12);
    constraintVar = [mu;gamma;phaseAng];
    constraintVar = constraintVar + J\M;
    mu = constraintVar(1:6);
    gamma = constraintVar(7:12);
    phaseAng = constraintVar(13);
    
    M
    pause
    flag2 = 0;
    for i = 1:13
        if abs(M(i)) > tolerance
            flag2 = 1;
            break;
        end
    end
%     
%     for i = 1:6
%         if mu(i) > pi/2
%             mu(i) = pi/3;
%         end
%     end
%     for i = 1:6
%         if gamma(i) > pi/2
%             gamma(i) = pi/3;
%         end
%     end
    
    if flag2 == 0
        break;
    end
    
    countIter = countIter + 1;
end

M
%pause

Pload = Vdc*Idc;

[Ia Ib Ic fundCurrents] = getI_abc_fundCurrents(I,h);