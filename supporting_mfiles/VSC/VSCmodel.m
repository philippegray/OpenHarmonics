%VSCmodel Returns the injection current harmonics for the converter of type
%"vsc".
function [Ia Ib Ic fundCurrents ma sigma Rdc_new error] = VSCmodel(Valphabeta,L,R,mf,Rdc,C,...
    Idc,h,Vdcref,Qload)
global OpenDSSFileLoc
Rdc_new = Rdc;
%L = 0.2;
% R = 0.02;
% Rdc = 10^8;
% C = 0.0181;
% Idc = 0.04;
% h = 15;
% mf = 15;
% Vdcref = 10;
% Qload = 1;
% Valphabeta = [0.0000;0.0000;0.0000;0.0000;0.0013;-0.0046;-0.0000;0.0000;0.0345;...
%         -0.0100;0.0000;0.0000;-0.0002;0.0011;0.0000;0.0000;-0.0244;...
%         -0.0336;0.0000;0.0000;-0.0295;-0.0811;0.0000;0.0000;0.0005;...
%         -0.0004;0.0000;0.0000;-0.0055;0.0027;0;0;0.9936;-0.0480;-0.0000;...
%         -0.0000;-0.0000;0.0003;0.0000;-0.0000;0.0152;-0.0133;-0.0000;...
%         0.0000;-0.0431;0.0651;0.0000;0.0000;0.0003;0.0011;-0.0000;0.0000;...
%         0.0050;-0.0019;-0.0000;-0.0000;0.0004;-0.0198;0.0000;0.0000;-0.0000;...
%         -0.0000];

w = 1; %system solved at 1 rad/s. System is per-unitized.
fac = w/2/pi;
m = 0; %not solving for any dc-harmonic terms
Pload = Vdcref*Idc; %0.2*2.73; %units of Watts

K = 3/2/(R^2+L^2);

% if abs(Pload) < 0.01
%     Idc = -0.001;
% end
if Pload ~= 0
	sigma = atan(Pload/(Qload-K*L));
    ma = -Pload/(K*L*sin(sigma)*Vdcref/2);
else
    sigma = 0;
    ma = 0.8;
end

% ma = 0.88; %modulation index value assumed for the PWM
% sigma = -5.1*pi/180;%phase angle offset value assumed for the PWM
wac = w;
tolerance = 0.001; %tolerance percent set to 1%
loopCount = 1;
Tac = 1/fac;
error = 0;
errorFile = fopen([OpenDSSFileLoc,'errorFile.txt'], 'w');
overModCount = 0;
%Now onto the calculation of PHI
% swStates = [0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,0;...
%     0,0,1,1,1,1,0,1,1,1,1,0,0,0,0,1,0,0,0;...
%     0,0,0,1,0,0,0,0,1,1,1,1,0,1,1,1,1,0,0];
% switchTimes = [0.0526,0.5303,0.9935,1.0998,1.5775,2.0407,2.1470,2.6247,3.0879,3.1942,3.6719,4.1351,4.2414,4.7191,5.1823,5.2886,5.7663,6.2295];
% numSWTimes = length(switchTimes);
%givens  **********

%Setting interface reactance of each phase to the input values.
Ra = R;
Rb = R;
Rc = R;
La = L;
Lb = L;
Lc = L;

%***********end givens
lenHarm = 2*h+1; %constant value for easy reference.

%This is the Clarke Transform. Used to convert between space-vector and abc
%reference frame.
CTf = 2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2)];
invCTf = inv(CTf);
Lmatrix = [La 0 0;0 Lb 0;0 0 Lc];
Rmatrix = [Ra 0 0;0 Rb 0;0 0 Rc];
%This is the symmetric matrix. Used to convert from abc to the
%positive-negative-zero sequence reference frame.
SymMatrix = [1 1 1;1 exp(1i*240*pi/180) exp(1i*120*pi/180);...
    1 exp(1i*120*pi/180) exp(1i*240*pi/180)];
invSymMatrix = inv(SymMatrix);

Mtx = transpose(SymMatrix)*conj(SymMatrix);
V = zeros(2*lenHarm+2*m+2,1);
V(2*lenHarm+1) = -Idc; %This is Idc, -Idc => that power is being drawn from grid.

V(1:2*lenHarm) = Valphabeta;
%voltageVec = V(1:2*lenHarm);
%save('Vvalues','voltageVec');

%The PWM scheme phase angle is with respect to the positive-sequence fundamental
%frequency component of the PCC voltages. So, phaseOff stores the phase
%angle of this fundamental frequecy component.
phaseOff = angle(V(lenHarm+2)+1i*V(lenHarm+3));

%Ts stores each time instant that each of the phases abc intersect with the
%control PWM sawtooth waveform. From these switching times we can
%reconstruct the state of the switches of the VSC at every instant through
%a whole period. From this we know the operation of the converter through a
%whole steady state period and can then fully model it.
Ts = zeros(2*mf+1,3); % Ts(1,:) is phase a, Ts(2,:) is phase b, Ts(3,:) is phase c
phaseArr = zeros(2*mf,3);
RotationMatrix = zeros(2*lenHarm+2*m+2,2*lenHarm+2*m+2);

while(true)
    mi = [ma;ma;ma];
    Za = 0+sigma+phaseOff; %phase offset of phase a at ti = 0 using the assumed sigma value. eac
    Zb = 4*pi/3+sigma+phaseOff; %phase offset of phase b at ti = 0 using the assumed sigma value.
    Zc = 2*pi/3+sigma+phaseOff; %phase offset of phase c at ti = 0 using the assumed sigma value.
    Zoff = [Za;Zb;Zc];
    
    %initializations for the loop
    %first period will be +1 to -1...then -1 to + 1...this continues...
    %the PWM switching times are calculated assuming the sample-and-hold
    %method.
    for k = 1:3 %this loop has three iterations-one iteration for each of the abc phases
        modphase = mi(k);
        for i = 1:2*mf+1%there will be 2*mf+1 switching times from ti -> ti + 2*\pi
            phase = Tac*(i-1)/(2*mf)+Zoff(k); %the point on 
            if mod(i,2) == 1
                swtime = fminbnd(@(x)pwmSwitchneg(x,modphase,phase,mf,Tac),0,Tac/2/mf);
                %swtime = Tac/4/mf*(1-ma*cos(phase));
            else
                swtime = fminbnd(@(x)pwmSwitchpos(x,modphase,phase,mf,Tac),0,Tac/2/mf);
                %swtime = Tac/4/mf*(1+ma*cos(phase));
            end
            Ts(i,k) = Tac/(2*mf)*(i-1)+swtime;
            phaseArr(i,k) = phase;
        end
    end
    
    data = [Ts(:,1),phaseArr(:,1),ones(2*mf+1,1);...
        Ts(:,2),phaseArr(:,2),2*ones(2*mf+1,1);...
        Ts(:,3),phaseArr(:,3),3*ones(2*mf+1,1)];
    data = sortrows(data);
    
    ti = data(1,1);
    numSWTimes = 6*mf;
    switchTimes = zeros(6*mf,1);
    for i = 1:6*mf
        switchTimes(i) = data(i+1,1)-data(i,1);
    end

    sw_data = zeros(6*mf,5);
    Sabc = zeros(3,1);
    for j = 1:mf
        Sabc = [0;0;0];
        for i = 0:1
            for k = 1:3
                posit = 6*(j-1)+3*i+k;
                Sabc(data(posit,3)) = mod(i+1,2);
                sw_data(posit,1:3) = transpose(Sabc);
            end
            Sabc = [1;1;1];
        end
    end
    for i = 1:6*mf
       sw_data(i,4:5) = [data(i,2),data(i+1,2)]; 
    end
    
    count = 1;
    for i  = -h:1:h
        RotationMatrix(count:count+1,count:count+1) = ...
            [cos(i*ti),-sin(i*ti);...
            sin(i*ti),cos(i*ti)];
        count = count + 2;
    end
    for i = 0:m
        RotationMatrix(count:count+1,count:count+1) = ...
            [cos(i*ti),-sin(i*ti);...
            sin(i*ti),cos(i*ti)];
        count = count + 2;
    end
    Zi = RotationMatrix*V;
    
    A11 = -CTf*inv(Lmatrix)*(Rmatrix*invCTf);
    A22 = -1/Rdc/C;
    
    A11 = A11(1:2,1:2);
    
    A = zeros(3,3,8);
    for Sc = 0:1
        for Sb = 0:1
            for Sa = 0:1
                index = Sc*2^2+Sb*2^1+Sa*2^0+1;
                A12 = -1/L*CTf(1:2,1:3)*[Sa;Sb;Sc];
                A21 = 1/C*[Sa Sb Sc]*invCTf(:,1:2);
                A(:,:,index) = [A11 A12;A21 A22];
            end
        end
    end
    
    Omega = zeros(2*lenHarm+2*m+2,2*lenHarm+2*m+2);
    M = zeros(2*lenHarm+2*m+2,2*lenHarm+2*m+2);
    %Calculation of N
    
    Nbase = CTf*inv(Lmatrix)*invCTf;
    Nbase = Nbase(1:2,1:2);
    N11 = zeros(2,lenHarm);
    N = zeros(3,2*lenHarm+2*m+2);
    
    for i = 0:lenHarm-1
        N11(1:2,2*i+1:2*i+2) = Nbase;
    end
    N(1:2,1:2*lenHarm) = N11;
    
    for i = 1:2:2*m+2
        N(3,2*lenHarm+i) = 1/C;
    end
    
    %Calculation of Omega
    sign = 1;
    count = 1;
    %Calculation of Omega = M
    for i = h:-1:0
        Omega(count:count+1,count:count+1) = [0 i*wac;-i*wac 0];
        count = count + 2;
    end
    for i = 1:h
        Omega(count:count+1,count:count+1) = [0 -i*wac;i*wac 0];
        count = count + 2;
    end
    for i = 0:m
        Omega(count:count+1,count:count+1) = [0 -i*wac;i*wac 0];
        count = count + 2;
    end
    
    M = Omega;
    H = zeros(2*lenHarm+2*m+2,3);
    
    %after re-deriving it looks lfine
    %Htemp = 1/Tac*exp(-1*1i*wac*Tac*harmArr'); %Re-add this later
    Htemp = 1/Tac*ones(lenHarm,1);
    
    count = 1;
    for i = h:-1:-h
        H(count:count+1,1:2) = 1/2/pi*[cos(i*2*pi) sin(i*2*pi);-sin(i*2*pi) cos(i*2*pi)];
        count = count + 2;
    end
    
    for i = 2*lenHarm+1:2:2*lenHarm+2*m+2
        H(i,3) = 1/Tac;
    end
    
    %Htemp = 3/pi*exp(-1*1i*pi/3*harmArr');
    % count = 1;
    % %Calculation of H
    % for i = 1:2*lenHarm
    %     if mod(i,2) == 0
    %         H(i,2) = Htemp(count);
    %         count = count+1;
    %     else
    %         H(i,1) = Htemp(count);
    %     end
    % end
    % Htemp = 0;%1/Tac*exp(1i*wac*Tac*harmArr'); %Re-add this later
    %
    % H(2*lenHarm+1:end,3) = Htemp;
    
    Ahat = zeros(3+2*lenHarm+2*m+2+2*lenHarm+2*m+2,3+2*lenHarm+2*m+2+2*lenHarm+2*m+2);
    Ahat(1:3,3+1:2*lenHarm+2*m+2+3) = N;
    Ahat(2*lenHarm+2*m+2+3+1:end,1:3) = H;
    Ahat(3+1:2*lenHarm+2*m+2+3,3+1:2*lenHarm+2*m+2+3) = Omega;
    Ahat(3+1+2*lenHarm+2*m+2:end,3+1+2*lenHarm+2*m+2:end) = M;
    %initializing phi*****
    phi = eye(length(Ahat),length(Ahat));
    leftFCMterms = zeros(length(phi),length(phi),numSWTimes-1);
    rightFCMterms = zeros(length(phi),length(phi),numSWTimes-1);
    
    index = sw_data(numSWTimes,3)*2^2+sw_data(numSWTimes,2)*2^1+sw_data(numSWTimes,1)*2^0+1;
    Ahat(1:3,1:3) = A(:,:,index);
    leftFCMterms(1:end,1:end,1) = expm(Ahat*switchTimes(end));
    
    index = sw_data(1,3)*2^2+sw_data(1,2)*2^1+sw_data(1,1)*2^0+1;
    Ahat(1:3,1:3) = A(:,:,index);
    rightFCMterms(1:end,1:end,end) = expm(Ahat*switchTimes(1));
    
    j = 2; %seems fine.
    for i = numSWTimes-1:-1:2 %might need to switch back to: 1:numSWTimes-1
        index = sw_data(i,3)*2^2+sw_data(i,2)*2^1+sw_data(i,1)*2^0+1;
        Ahat(1:3,1:3) = A(:,:,index);
        leftFCMterms(:,:,j) = leftFCMterms(:,:,j-1)*expm(Ahat*switchTimes(i));
        j = j + 1;
    end
    
    j = numSWTimes-2;
    for i = 2:numSWTimes-1
        index = sw_data(i,3)*2^2+sw_data(i,2)*2^1+sw_data(i,1)*2^0+1;
        Ahat(1:3,1:3) = A(:,:,index);
        rightFCMterms(:,:,j) = expm(Ahat*switchTimes(i))*rightFCMterms(:,:,j+1);
        j = j - 1;
    end
    
    phi = leftFCMterms(:,:,5)*rightFCMterms(:,:,5);
    
    %IF THE LAST SWITCH TIME IS LESS THEN the period Ts THEN NEED TO ADD
    %ANOTHER TERM BECAUSE ALL SWITCHES WILL BE ON until time = 0;
    
    Hp = phi(2*lenHarm+2*m+2+3+1:end,1:3);
    Ap = phi(1:3,1:3);
    Np = phi(1:3,3+1:2*lenHarm+2*m+2+3);
    Qp = phi(2*lenHarm+2*m+2+3+1:end,3+1:3+2*lenHarm+2*m+2);
    
    FCM = Hp*(inv(eye(3,3)-Ap))*Np+Qp;
    VIarray_alphabeta = zeros(h,4);
    
    %Testing Purposes*********************************************
    Ialphabeta = FCM*Zi;
    
    Vdc = Ialphabeta(2*lenHarm+1);
    
%     dPhi_dma = Tac/4/mf*(-cos(switchTimesData(1,2)))*...
%         leftFCMterms(1:end,1:end,numSWTimes)*Ahat*...
%         rightFCMterms(1:end,1:end,numSWTimes);
%     dPhi_dma = dPhi_dma + Tac/4/mf*Ahat*(cos(switchTimesData(numSWTimes,2)))*phi;
%     
%     dPhi_dsigma = Tac/4/mf*ma*(sin(switchTimesData(1,2)))*...
%         leftFCMterms(1:end,1:end,numSWTimes)*Ahat*...
%         rightFCMterms(1:end,1:end,numSWTimes);
%     dPhi_dsigma = dPhi_dsigma - Tac/4/mf*ma*Ahat*(sin(switchTimesData(numSWTimes,2)))*phi;    
    
    dPhi_dma = 0;
    dPhi_dsigma = 0;
    sub = 0;
    
    for i = 1:mf
        if i == mf
           sub = 1; 
        end
        for j = 1:6-sub
            pos = (i-1)*6+j;
            Ang1 = sw_data(pos,4);
            Ang2 = sw_data(pos,5);
            if j <= 2
                dmaTerm = Tac/4/mf*(cos(Ang1)-cos(Ang2));
                dsigmaTerm = Tac/4/mf*ma*(-sin(Ang1)+sin(Ang2));
            elseif j == 3
                dmaTerm = Tac/4/mf*(cos(Ang2)+cos(Ang1));
                dsigmaTerm = -Tac/4/mf*ma*(sin(Ang2)+sin(Ang1));
            elseif j <= 5
                dmaTerm = -Tac/4/mf*(cos(Ang1)-cos(Ang2));
                dsigmaTerm = -Tac/4/mf*ma*(-sin(Ang1)+sin(Ang2));
            else
                dmaTerm = -Tac/4/mf*(cos(Ang2)+cos(Ang1));
                dsigmaTerm = Tac/4/mf*ma*(sin(Ang2)+sin(Ang1));
            end
            index = sw_data(pos,3)*2^2+sw_data(pos,2)*2^1+sw_data(pos,1)*2^0+1;
            Ahat(1:3,1:3) = A(:,:,index);
            
            dPhi_dma = dPhi_dma + leftFCMterms(:,:,numSWTimes-pos)*...
                dmaTerm*Ahat*rightFCMterms(:,:,numSWTimes-pos);
            
            dPhi_dsigma = dPhi_dsigma + leftFCMterms(:,:,numSWTimes-pos)*...
                dsigmaTerm*Ahat*rightFCMterms(:,:,numSWTimes-pos);
        end
    end
    
    dHp_dma = dPhi_dma(2*lenHarm+2*m+2+3+1:end,1:3);
    dAp_dma = dPhi_dma(1:3,1:3);
    dNp_dma = dPhi_dma(1:3,3+1:2*lenHarm+2*m+2+3);
    dQp_dma = dPhi_dma(2*lenHarm+2*m+2+3+1:end,3+1:3+2*lenHarm+2*m+2);
    
    dHp_dsigma = dPhi_dsigma(2*lenHarm+2*m+2+3+1:end,1:3);
    dAp_dsigma = dPhi_dsigma(1:3,1:3);
    dNp_dsigma = dPhi_dsigma(1:3,3+1:2*lenHarm+2*m+2+3);
    dQp_dsigma = dPhi_dsigma(2*lenHarm+2*m+2+3+1:end,3+1:3+2*lenHarm+2*m+2);
    
    %THIS IS INCORRECT LOOK AT WHAT I DO FOR THE THYRISTOR
    %MODEL!!!!!!!!!!!!!
    dFCM_dma = derivFCM(dHp_dma,dNp_dma,dQp_dma,dAp_dma,Hp,Np,Ap);
    dFCM_dsigma = derivFCM(dHp_dsigma,dNp_dsigma,dQp_dsigma,dAp_dsigma,Hp,Np,Ap);%dHp,dNp,dQp,dAp,Hp,Np,Qp,Ap
    
    %**********************************************************************
    
    dQ_dma = 0;
    dQ_dsigma = 0;

    %To calculate dQ_dma and dQ_dsigma***
%     for i = 1:2:2*lenHarm
%         dIsa_dma = dFCM_dma(i,1:end)*Zi;
%         dIsb_dma = dFCM_dma(i+1,1:end)*Zi;
%         dIsa_dsigma = dFCM_dsigma(i,1:end)*Zi;
%         dIsb_dsigma = dFCM_dsigma(i+1,1:end)*Zi;
%         dQ_dma = dQ_dma-3/2*Zi(i)*dIsb_dma+3/2*Zi(i+1)*dIsa_dma;
%         dQ_dsigma = dQ_dsigma-3/2*Zi(i)*dIsb_dsigma+3/2*Zi(i+1)*dIsa_dsigma;
% %         dP_dma = dP_dma+3/2*V(i)*dIsa_dma+3/2*V(i+1)*dIsb_dma;
% %         dP_dsigma = dP_dsigma+3/2*V(i)*dIsa_dsigma+3/2*V(i+1)*dIsb_dsigma;
%     end
    i = 2*h+3;
    dIsa_dma = dFCM_dma(i,1:end)*Zi;
    dIsb_dma = dFCM_dma(i+1,1:end)*Zi;
    dIsa_dsigma = dFCM_dsigma(i,1:end)*Zi;
    dIsb_dsigma = dFCM_dsigma(i+1,1:end)*Zi;
    dQ_dma = dQ_dma-3/2*Zi(i)*dIsb_dma+3/2*Zi(i+1)*dIsa_dma;
    dQ_dsigma = dQ_dsigma-3/2*Zi(i)*dIsb_dsigma+3/2*Zi(i+1)*dIsa_dsigma;
    i = 2*h;
    dIsa_dma = dFCM_dma(i,1:end)*Zi;
    dIsb_dma = dFCM_dma(i+1,1:end)*Zi;
    dIsa_dsigma = dFCM_dsigma(i,1:end)*Zi;
    dIsb_dsigma = dFCM_dsigma(i+1,1:end)*Zi;
    dQ_dma = dQ_dma-3/2*Zi(i)*dIsb_dma+3/2*Zi(i+1)*dIsa_dma;
    dQ_dsigma = dQ_dsigma-3/2*Zi(i)*dIsb_dsigma+3/2*Zi(i+1)*dIsa_dsigma;    
    
    
    dP_dma = Idc*dFCM_dma(2*lenHarm+1,1:end)*Zi;
    dP_dsigma = Idc*dFCM_dsigma(2*lenHarm+1,1:end)*Zi;
    
%   Scalc = 0;
%     for i = 0:lenHarm-1
%         Scalc = Scalc+3/2*(Zi(2*i+1)+1i*Zi(2*i+2))*...
%             conj(Ialphabeta(2*i+1)+1i*Ialphabeta(2*i+2));
%   end
    Qnegseq = 3/2*imag((Zi(2*h+3)+1i*Zi(2*h+4))*conj(Ialphabeta(2*h+3)+1i*Ialphabeta(2*h+4)));
    Qposseq = 3/2*imag((Zi(2*h-1)+1i*Zi(2*h))*conj(Ialphabeta(2*h-1)+1i*Ialphabeta(2*h)));
    Qcalc = Qposseq+Qnegseq;

    %Qcalc = imag(Scalc);
    Pcalc = Idc*Vdc;
    Jload = [dP_dsigma,dP_dma;dQ_dsigma,dQ_dma];
    
    delP = Pload - Pcalc;
    delQ = Qload - Qcalc;
    
    temp = [sigma;ma] + Jload\[delP;delQ];
    sigma = temp(1);
    ma = temp(2);
    
    if (ma > 1)
        overModCount = overModCount + 1;
        if overModCount > 4
            c = clock;
            fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Solution does not exist. Increase Vdc.\n',c(4),c(5),c(6));
            error = 1;
        	break;
        end
    elseif overModCount > 0
       overModCount = 0; 
    end
    
%     if (ma > 2)
%         c = clock;
%         fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Solution does not exist. Increase Vdc.\n',c(4),c(5),c(6));
%         error = 1;
%     end
    
    if ma < 0
        ma = 0.01;
    end
%     if (sigma > 2*pi)
%         sigma = 0;
%     end
%     if sigma < -2*pi
%         sigma = 0;
%     end
    
    if abs(delP) < tolerance && abs(delQ) < tolerance
        break
    end
    loopCount = loopCount+1;
    
    if loopCount == 50
        c = clock;
        fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Solution does not exist for VSC\n',c(4),c(5),c(6));
        error = 1;
        break;
    end
end

if ma > 1 && overModCount <= 4 
    c = clock;
    fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Solution does not exist. Increase Vdc.\n',c(4),c(5),c(6));
    error = 1;
end
        
count = 1;
for i  = -h:1:h
    RotationMatrix(count:count+1,count:count+1) = ...
        [cos(-i*ti),-sin(-i*ti);sin(-i*ti) cos(-i*ti)];
    count = count + 2;
end
for i  = 0:1:m
    RotationMatrix(count:count+1,count:count+1) = ...
        [cos(-i*ti),-sin(-i*ti);sin(-i*ti) cos(-i*ti)];
    count = count + 2;
end

%inv(RotationMatrix)*
tempM = eye(2*lenHarm+2*m+2);
tempM(2*lenHarm+3:end,2*lenHarm+1:end) = 2*tempM(2*lenHarm+3:end,2*lenHarm+1:end);
Ialphabeta = tempM*RotationMatrix*Ialphabeta;
sigma = mod(sigma+phaseOff,2*pi);
[Ia Ib Ic fundCurrents] = getI_abc_fundCurrents(Ialphabeta,h);
%plotFCMaxes(visualizeFCM,h);
%plotRXaxes(iFCMRXmag,h);
%plotMNaxes(iFCMMNmag,h);

fclose(errorFile);
