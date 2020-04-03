function [Ia Ib Ic fundCurrents ma sigma] = VSCmodel(Valphabeta,L,R,mf,Rdc,C,...
    Idc,h,Vdcref,Qload)

%global SymMtx

%function 
%Givens
% fac = 60;
% T = 1/fac;
% wac = 2*pi*fac;

w = 1;
fac = w/2/pi;
m = 0;
Pload = Vdcref*Idc;%0.2*2.73; %units of Watts
ma = 0.88; %modulation index of a
sigma = -5.1*pi/180;%-5*pi/180;
wac = w;
tolerance = 1;
loopCount = 1;

%sigma = -8*pi/180;%-1.2*pi/180;%pi/6;%Initially setting sigma to 0...will iterate to find value of
%sigma that gives output power = Pload%-1.2*pi/180;
Tac = 1/fac;

sizeArr = 10000;

%Now onto the calculation of PHI
% swStates = [0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,0;...
%     0,0,1,1,1,1,0,1,1,1,1,0,0,0,0,1,0,0,0;...
%     0,0,0,1,0,0,0,0,1,1,1,1,0,1,1,1,1,0,0];
% switchTimes = [0.0526,0.5303,0.9935,1.0998,1.5775,2.0407,2.1470,2.6247,3.0879,3.1942,3.6719,4.1351,4.2414,4.7191,5.1823,5.2886,5.7663,6.2295];
% numSWTimes = length(switchTimes);
%givens  **********
Ra = R;
Rb = R;
Rc = R;
La = L;
Lb = L;
Lc = L;
w = 1;

%***********end givens
harmArr = -h:1:h;
lenHarm = 2*h+1;

%Rmatrix = [Raa Rab Rac;Rba Rbb Rbc;Rca Rcb Rcc];

CTf = 2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2)];
invCTf = inv(CTf);
Lmatrix = [La 0 0;0 Lb 0;0 0 Lc];
Rmatrix = [Ra 0 0;0 Rb 0;0 0 Rc];
L = La;
R = Ra;
%xabc = SymMatrix*x0+-
SymMatrix = [1 1 1;1 exp(1i*240*pi/180) exp(1i*120*pi/180);...
    1 exp(1i*120*pi/180) exp(1i*240*pi/180)];
invSymMatrix = inv(SymMatrix);
X = w*L;

V = zeros(2*lenHarm+2*m+2,1);
V(2*lenHarm+1) = -Idc; %This is Idc, -Idc => that power is being drawn from grid.

V(1:2*lenHarm) = Valphabeta;
%voltageVec = V(1:2*lenHarm);
%save('Vvalues','voltageVec');

phaseOff = angle(V(lenHarm+2)+1i*V(lenHarm+3));
%Initialization

Ts = zeros(3,mf); % Ts(1,:) is phase a, Ts(2,:) is phase b, Ts(3,:) is phase c

while(true)
    mi = [ma;ma;ma];
    Za = 0+sigma+phaseOff; %phase offset of phase a
    Zb = 4*pi/3+sigma+phaseOff; %phase offset of phase b
    Zc = 2*pi/3+sigma+phaseOff; %phase offset of phase c
    Zoff = [Za;Zb;Zc];
    
    %initializations for the loop
    %first period will be +1 to -1...then -1 to + 1...this continues...
    for k = 1:3
        modphase = mi(k);
        for i = 1:2*mf
            phase = Tac*(i-1)/(2*mf)+Zoff(k);
            if mod(i,2) == 1
                swtime = fminbnd(@(x)pwmSwitchneg(x,modphase,phase,mf,Tac),0,Tac/2/mf);
            else
                swtime = fminbnd(@(x)pwmSwitchpos(x,modphase,phase,mf,Tac),0,Tac/2/mf);
            end
            Ts(k,i) = Tac/(2*mf)*(i-1)+swtime;
        end
    end
    
    switchTimes = sort([Ts(1,:),Ts(2,:),Ts(3,:)]); %switching times are correct for Sample and Hold method!
    
    numSWTimes = length(switchTimes);
    swStates = zeros(3,numSWTimes+1);
    
    %For determining at which times the switches of the various phases are on
    %or off.
    for k = 1:3 %something is wrong here.
        count = 1;
        currState = -1;
        
        for i = 1:numSWTimes
            if (count <= mf*2 && switchTimes(i) == Ts(k,count))
                currState = currState*-1;
                count = count+1;
            end
            
            if (currState == -1)
                swStates(k,i+1) = 0;
            else
                swStates(k,i+1) = currState;
            end
        end
    end
    
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
                %A21 = [A21(1),A21(2)];
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
    Ahat(1:3,1:3) = A(:,:,1);
    phi = expm(Ahat*(switchTimes(1)));
    %*****
    
    for i = 2:(numSWTimes) %might need to switch back to: 1:numSWTimes-1
        index = swStates(3,i)*2^2+swStates(2,i)*2^1+swStates(1,i)*2^0+1;
        Ahat(1:3,1:3) = A(:,:,index);
        phi = expm(Ahat*(switchTimes(i)-switchTimes(i-1)))*phi;
    end
    %IF THE LAST SWITCH TIME IS LESS THEN the period Ts THEN NEED TO ADD
    %ANOTHER TERM BECAUSE ALL SWITCHES WILL BE ON until time = 0;
    
    Ahat(1:3,1:3) = A(:,:,1);
    phi = expm(Ahat*(Tac-switchTimes(numSWTimes)))*phi;
    
    Hp = phi(2*lenHarm+2*m+2+3+1:end,1:3);
    Ap = phi(1:3,1:3);
    Np = phi(1:3,3+1:2*lenHarm+2*m+2+3);
    Qp = phi(2*lenHarm+2*m+2+3+1:end,3+1:3+2*lenHarm+2*m+2);
    
    %FCM = abs(Hp*(inv(eye(3,3)-Ap))*Np+Qp);
    FCM = Hp*(inv(eye(3,3)-Ap))*Np+Qp;
    VIarray_alphabeta = zeros(h,4);
    
    %Testing Purposes*********************************************
    Ialphabeta = FCM*V;
    
    Vval = abs(V(lenHarm+2)+1i*V(lenHarm+3));
    Vdc = Ialphabeta(2*lenHarm+1);
    
    dQ_dSigma = 3/2*Vval/(R^2+X^2)*(ma*Vdc/2*R*cos(sigma)+ma*Vdc/2*X*sin(sigma));
    dQ_dma = 3/2*Vval/(R^2+X^2)*(Vdc/2*R*sin(sigma)-Vdc/2*X*cos(sigma));
    %dV_dSigma = -2/ma*Vval*sin(sigma)-2/ma*R*X*Vval*cos(sigma);
    %dV_dma = -2/ma^2*Vval*cos(sigma)+2/ma^2*R*X*Vval*sin(sigma)+16/3/ma^3/R*Idc*(R^2+X^2);
    %dQ_dSigma = w*C*2*Vdc*dV_dSigma;
    %dQ_dma = w*C*2*Vdc*dV_dma;
    dP_dSigma = 3/2*Vval/(R^2+X^2)*(ma*Vdc/2*R*sin(sigma)-ma*Vdc/2*X*cos(sigma));
    dP_dma = -3/2*Vval/(R^2+X^2)*(Vdc/2*R*cos(sigma)+Vdc/2*X*sin(sigma));
    %Pcalc = Vdc*Idc;%3/4/(R^2+X^2)*Vdc*ma*(Vval*cos(sigma)*R-Vval*X*sin(sigma)-Vdc/2*ma*R);
    %Pcalc =
    %3*Vdc*ma/2/(R^2+X^2)*(Vval*R*cos(sigma)-Vdc/2*ma*R-Vval*X*sin(sigma));
    %Qcalc = 3*Vdc*ma/2/(R^2+X^2)*(Vval*X*cos(sigma)-X*Vdc/2*ma+Vval*R*sin(sigma));
    
    %Testing****
    Scalc = 0;
    for i = 0:lenHarm-1
        Scalc = Scalc+3/2*(V(2*i+1)+1i*V(2*i+2))*...
            conj(Ialphabeta(2*i+1)+1i*Ialphabeta(2*i+2));
    end
    %Pcalc = real((V(lenHarm+2)+1i*V(lenHarm+3))*conj(Ialphabeta(lenHarm+2)+1i*Ialphabeta(lenHarm+3)))
    %-3*R*abs(Ialphabeta(lenHarm+2)+1i*Ialphabeta(lenHarm+3))^2;
    
    %***********
    Qcalc = imag(Scalc);
    Pcalc = real(Scalc);
    Pcalc = Idc*Vdc;
    Jload = [dP_dSigma,dP_dma;dQ_dSigma,dQ_dma];
    
    delP = Pload - Pcalc;
    delQ = Qload - Qcalc;
    
    temp = [sigma;ma] + Jload\[delP;delQ];
    sigma = temp(1);
    ma = temp(2);
    
    if (ma > 1)
        ma = 0.5;
    end
    if ma < 0
        ma = 0.5;
    end
    
    if sqrt((delP/Pload*100)^2+(delQ/Qload*100)^2) < tolerance
        break
    end
    loopCount = loopCount+1;
    
    if loopCount == 50
        display('Solution cannot be found! Specify a new PQ value!');
        break
    end
end
loopCount
tempM = eye(2*lenHarm+2*m+2);
tempM(2*lenHarm+3:end,2*lenHarm+1:end) = 2*tempM(2*lenHarm+3:end,2*lenHarm+1:end);
Ialphabeta = tempM*Ialphabeta;

[Ia Ib Ic fundCurrents] = getI_abc_fundCurrents(Ialphabeta,h);
%plotFCMaxes(visualizeFCM,h);
%plotRXaxes(iFCMRXmag,h);
%plotMNaxes(iFCMMNmag,h);
