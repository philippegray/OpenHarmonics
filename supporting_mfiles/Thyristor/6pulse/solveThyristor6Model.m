function [Ia,Ib,Ic,fundCurrents,alpha,mu_out,error] = solveThyristor6Model(Valphabeta,L,R,...
    Rdc,Ldc,h,Vdc,P,mu_in,alpha_in)

global SymMtx invSymMtx


countIter = 0;
Pload = P;
Rout = Rdc;
Lout = Ldc;
w = 1;     %line frequency in rad/s
alpha = alpha_in;
tolerance = 0.001;
f = w/2/pi;
Tac = 1/f;
Ra = R;
Rb = R;
Rc = R;
La = L;
Lb = L;
Lc = L;
lenHarm = 2*h+1;
%Mag_abc = [0.8; 1;0.8]; %units: V_{pk}
%Phase_abc = [0/180*pi+0*pi/180;240*pi/180+0*pi/180;120*pi/180+0*pi/180]; %why was every phase angle before phase shifted by -90 rad.
%V_abc = A*[0.1;1+1*1i;0.05+0.05*1i];
V_abc = SymMtx*[0;Valphabeta(lenHarm+2)+1i*Valphabeta(lenHarm+3);...
    Valphabeta(2*h-1)+1i*Valphabeta(2*h)];
Mag_abc = [abs(V_abc(1));abs(V_abc(2));abs(V_abc(3))];
Phase_abc = [angle(V_abc(1));angle(V_abc(2));angle(V_abc(3))];
error = 0;
CTF = 2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2)];
invCTF = inv(CTF);

Mi = zeros(6,2);
Mi(2,1:end) = [1 0];
Mi(3,1:end) = -1*[-1/2 -sqrt(3)/2];
Mi(1,1:end) = -1*[-1/2 sqrt(3)/2];
Mi(5,1:end) = -1*[1 0];
Mi(6,1:end) = [-1/2 -sqrt(3)/2];
Mi(4,1:end) = [-1/2 sqrt(3)/2];
%end initializations************************

Vsabc_phasor = zeros(3,1);
for i = 1:3
    [x y] = pol2cart(Phase_abc(i),Mag_abc(i));
    Vsabc_phasor(i) = x + 1i*y;
end

temp = [1 -1 0;0 1 -1;-1 0 1]*[Vsabc_phasor(1);Vsabc_phasor(2);Vsabc_phasor(3)];
    
Mag_sv = abs(temp);
Phase_sv = angle(temp);
flag_error = 0;

maxPos = 3+3*lenHarm+1;

%Testing out new algorithm for switching times**********************
symMatrix = invSymMtx*Vsabc_phasor;

posSeq = symMatrix(2);
negSeq = symMatrix(3);

phase = angle(posSeq)-0;%-30/180*pi;%150/180*pi;
Ts = zeros(18,1);
Ts(1:end) = [1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6];

Zi = zeros(3*lenHarm+1,1);
Zi(1:2*lenHarm,1) = Valphabeta;
Zi(2*lenHarm+1) = Vdc; %This is Idc, -Idc => that power is being drawn from grid.
%**************************************************************************

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

mu = zeros(7,1); %Initialization
mu(1:6) = mu_in;
mu(7) = alpha_in;
%The following for loop calculates the delta mu(1)...mu(3) terms in the
%Jacobian matrix formulation.
new_mu = zeros(7,1);
tempArr = [1 0 0;0 1 0;0 0 1];
At = zeros(3,3);
Nt = zeros(3,4*lenHarm);
Omegat = zeros(3*lenHarm+1,3*lenHarm+1);
ThetaM = [1 0;0 1];

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

Mcond = zeros(3+3*lenHarm+1,3+3*lenHarm+1);
Mcomm = zeros(3+3*lenHarm+1,3+3*lenHarm+1);
Mt_arr = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6);
derivMp = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6,6); %3rd column contains the partial
Mp = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6);
Ap = zeros(3,3,6);
Np = zeros(3,3*lenHarm+1,6);
Omegap = zeros(3*lenHarm+1,3*lenHarm+1,6);

%This loop is used to solve for the alpha that will yield the required
%power sourced by the load.

Ipwr_alpha = 0;
Ipwr_beta = 0;

Ahat = zeros(3+3*lenHarm+1+3*lenHarm+1);

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

Ahat(3+3*lenHarm+1+1:end,1:3) = H;
Ahat(4:3+3*lenHarm+1,4:3+3*lenHarm+1) = Omegat;
Ahat(3+3*lenHarm+1+1:end,3+3*lenHarm+1+1:end) = Mm;

while(true)
    
    
    %This loop is used to solve for mu(1) --> mu(6) through a Newton Raphson
    %Iterative solving technique.
    
    RotationMatrix = zeros(3*lenHarm+1,3*lenHarm+1,6);
    Zi_arr = zeros(3*lenHarm+1,6);
    for k = 0:1:5
        count = 1;
        for i  = -h:1:h
            RotationMatrix(count:count+1,count:count+1,k+1) = ...
                [cos(i*(k*pi/3+alpha-phase)),...
                -sin(i*(k*pi/3+alpha-phase));...
                sin(i*(k*pi/3+alpha-phase)),...
                cos(i*(k*pi/3+alpha-phase))];
            count = count + 2;
        end
        for i = 0:h
            RotationMatrix(count:count+1,count:count+1,k+1) = ...
                [cos(i*(k*pi/3+alpha-phase)),...
                -sin(i*(k*pi/3+alpha-phase));...
                sin(i*(k*pi/3+alpha-phase)),...
                cos(i*(k*pi/3+alpha-phase))];
            count = count + 2;
        end
        Zi_arr(1:end,k+1) = RotationMatrix(1:end,1:end,k+1)*Zi;
    end
    
    %This block of code generates the Mp matrix elements and the dMp/dx
    %elements
    
    %This block of code obtaines all the matrix exponentials that will be
    %needed for the first part of the code.
    
    %third column is referring to interval with respect to the primary
    %rotation configuration
    %fourth column contains 3 particular exponential calculations
    expmCommMu = zeros(length(Ahat),length(Ahat),6);
    expmCondMu = zeros(length(Ahat),length(Ahat),6);
    expmCondMuPi_3 = zeros(length(Ahat),length(Ahat),6);
    
    Phi = eye(3+3*lenHarm+1+3*lenHarm+1);
    for i = 1:6
        Ahat(1:3,1:3+3*lenHarm+1) = [Acomm(1:end,1:end,i),Ncomm(1:end,1:end,i)];
        expmCommMu(1:end,1:end,i) = expm(Ahat*mu(i));
        Ahat(1:3,1:3+3*lenHarm+1) = [Acond(1:end,1:end,i),Ncond(1:end,1:end,i)];
        expmCondMu(1:end,1:end,i) = expm(-Ahat*mu(i));
        expmCondMuPi_3(1:end,1:end,i) = expm(Ahat*pi/3);       
    
        Phi = expmCondMuPi_3(1:end,1:end,i)*...
            expmCondMu(1:end,1:end,i)*expmCommMu(1:end,1:end,i)*Phi;
    end
    
    %Calculating FCM*******************************************************
    Hp = Phi(4+3*lenHarm+1:end,1:3);
    Ap_mat = Phi(1:3,1:3);
    Np_mat = Phi(1:3,4:3*lenHarm+4);
    Qp = Phi(3*lenHarm+4+1:end,4:3*lenHarm+4);
    
    FCM = Hp*(inv(eye(3,3)-Ap_mat))*Np_mat+Qp;
    I = FCM*Zi_arr(1:end,1);
    %**********************************************************************
    for i = 1:6
        Mp(1:end,1:end,i) = eye(length(Mp));
        for j = 1:6
            derivMp(1:end,1:end,i,j) = eye(length(derivMp));
        end
    end
    %Calculating Mp and derivMp********************************************
    for i = 1:6 %i refers to M1 through to M6
        for j = 0:5 %j steps through to account for a whole period
            Mcond(1:3,1:3) = Acond(1:end,1:end,Ts(i+j));
            Mcond(1:3,4:end) = Ncond(1:end,1:end,Ts(i+j));
            Mcond(4:end,4:end) = Omegat;
            Mcomm(1:3,1:3) = Acomm(1:end,1:end,Ts(i+j));
            Mcomm(1:3,4:end) = Ncomm(1:end,1:end,Ts(i+j));
            Mcomm(4:end,4:end) = Omegat;
            
            %Mp(*,*,1) goes from ti-phase+alpha -> ti-phase+alpha+2*pi
            %Mp(*,*,3) goes from ti-phase+alpha+2pi/3 ->
            %ti-phase+alpha+8pi/3
            Mp(1:end,1:end,i) = expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*Mp(1:end,1:end,i);
            
            for k = 1:6 %k refers to taking derivative wrt commutation overlap length with indice k
                if k == (j+1)
                    derivMp(1:end,1:end,k,i) = expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                        (-Mcond*expmCondMu(1:maxPos,1:maxPos,Ts(i+j))+...
                        expmCondMu(1:maxPos,1:maxPos,Ts(i+j))*Mcomm)*expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                        derivMp(1:end,1:end,k,i);
                else
                    derivMp(1:end,1:end,k,i) = expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                        expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*derivMp(1:end,1:end,k,i);
                end
            end
        end
    end
    
    dMp_dmu = zeros(length(Mp),length(Mp),6,6);
    
    for i = 0:1:5
        for j = 1:1:6
            dMp_dmu(1:end,1:end,Ts(j+i),i+1) = derivMp(1:end,1:end,j,i+1);
        end
    end
    %derivMp = dMp_dmu;
    %**********************************************************************
    
    %the 4th dimension is whether it be M1,M2,M3 (each has an offset to the initial vector)
    %the 3rd dimension is for the derivative with respect to mu1,mu2,mu3
    
    %[x;u]@t+T = Mp(*,*,1)*[x;u]@t
    for i = 1:6
        Ap(1:end,1:end,i) = Mp(1:3,1:3,i);
        Np(1:end,1:end,i) = Mp(1:3,4:end,i);
        Omegap(1:end,1:end,i) = Mp(4:end,4:end,i);
    end
    
    derivMt = zeros(3+3*lenHarm+1,3+3*lenHarm+1,6);
    
    %M2 SOLVE
    %This block calculates the derivative of the Mt matrix elements.*******
    for i = 1:6
        Mcomm(1:3,1:3) = Acomm(1:end,1:end,i);
        Mcomm(1:3,4:end) = Ncomm(1:end,1:end,i);
        Mcomm(4:end,4:end) = Omegat;
        
        Mt_arr(1:end,1:end,i) = expmCommMu(1:maxPos,1:maxPos,i);
        derivMt(1:end,1:end,i) = Mcomm*Mt_arr(1:end,1:end,i);
    end
    %**********************************************************************
    
    dPhi_dmu = zeros(length(Phi),length(Phi),6);
    Ahat2 = Ahat;
    
    for i = 1:6
        dPhi_dmu(1:end,1:end,i) = eye(length(Phi));
        for j = 1:6
            if i == j
                Ahat(1:3,1:3+3*lenHarm+1) = [Acomm(1:end,1:end,j),Ncomm(1:end,1:end,j)];
                Ahat2(1:3,1:3+3*lenHarm+1) = [Acond(1:end,1:end,j),Ncond(1:end,1:end,j)];
                dPhi_dmu(1:end,1:end,i) = expmCondMuPi_3(1:end,1:end,j)*...
                    (-Ahat2*expmCondMu(1:end,1:end,j)+expmCondMu(1:end,1:end,j)*Ahat)...
                    *expmCommMu(1:end,1:end,j)*dPhi_dmu(1:end,1:end,i);
            else
                dPhi_dmu(1:end,1:end,i) = expmCondMuPi_3(1:end,1:end,j)*...
                    expmCondMu(1:end,1:end,j)*expmCommMu(1:end,1:end,j)*...
                    dPhi_dmu(1:end,1:end,i);
            end
        end
    end
    
    derivRotation = zeros(length(FCM),length(FCM),6);
    for j = 0:1:5
        count = 1;
        for i  = -h:1:h
            derivRotation(count:count+1,count:count+1,j+1) = [-i*sin(i*(j*pi/3+alpha-phase)),...
                -i*cos(i*(j*pi/3+alpha-phase));i*cos(i*(j*pi/3+alpha-phase)) -i*sin(i*(j*pi/3+alpha-phase))];
            count = count + 2;
        end
        for i  = 0:1:h
            derivRotation(count:count+1,count:count+1,j+1) = [-i*sin(i*(j*pi/3+alpha-phase)),...
                -i*cos(i*(j*pi/3+alpha-phase));i*cos(i*(j*pi/3+alpha-phase)) -i*sin(i*(j*pi/3+alpha-phase))];
            count = count + 2;
        end
    end
    
    Idc = I(2*lenHarm+1);
    Vdc = Zi(2*lenHarm+1);
    Irow_FCM_dalpha = FCM*derivRotation(1:end,1:end,1);
    Irow_FCM_dalpha = Irow_FCM_dalpha(2*lenHarm+1,1:end);
    dIdc_dalpha = Irow_FCM_dalpha*Zi;
    
    %I also need the correct Mcomm...
    %There has to be a conditional if statement that if j ==
    %arrPosArr(Ts(6,2)), then for M1 M2 M3 the values have to be saved for
    %that derivative term.  Actually, might not even need to use a loop for
    %this...could just directly add the terms as no further derivateves are
    %required.
    %**************************************************
    J = zeros(7,7);
    for i = 1:6
        At = Mt_arr(1:2,1:2,i);
        Nt = Mt_arr(1:2,4:end,i);
        J(i,7) = Mi(i,1:end)*(At*inv(eye(2)-Ap(1:2,1:2,i))*...
            Np(1:2,1:end,i)+Nt)*derivRotation(1:end,1:end,i)*Zi;
        for j = 1:6
            if i == j
                dAt = derivMt(1:2,1:2,i);
                dNt = derivMt(1:2,4:end,i);
            else
                dAt = zeros(2,2);
                dNt = zeros(2,3*lenHarm+1);
            end
            dAp = dMp_dmu(1:2,1:2,j,i);
            dNp = dMp_dmu(1:2,4:end,j,i);
            invA = inv(eye(2)-Ap(1:2,1:2,i));
            det = (1 - Ap(1,1,i))*(1-Ap(2,2,i))-Ap(1,2,i)*Ap(2,1,i);
            ddet_dmu = (-1+Ap(2,2,i))*dAp(1,1)-Ap(2,1,i)*dAp(1,2)-...
                Ap(1,2,i)*dAp(2,1)+(Ap(1,1)-1)*dAp(2,2);
            derivInv = 1/(det^2)*[-det*dAp(2,2)-(1-Ap(2,2,i))*ddet_dmu,...
               det*dAp(1,2)-Ap(1,2,i)*ddet_dmu;det*dAp(2,1)-Ap(2,1,i)*ddet_dmu,...
               -det*dAp(1,1)-(1-Ap(1,1,i))*ddet_dmu];
            J(i,j) = Mi(i,1:end)*(dAt*invA*Np(1:2,1:end,i)+At*invA*dNp+...
                At*derivInv*Np(1:2,1:end,i)+dNt)*Zi_arr(1:end,i);
        end
    end
         
    dFCM_dmu = zeros(length(FCM),length(FCM),6);
    derivInv2 = zeros(3,3);
    derivAp = zeros(3,3);
    
    %Working now
    for i = 1:6 %FCM = Hp*(inv(eye(3,3)-Ap_mat))*Np_mat+Qp;
        dAp = dPhi_dmu(1:3,1:3,i);
        
        det = (Ap_mat(3,3)-1)*(Ap_mat(1,1)+Ap_mat(2,2)-Ap_mat(1,1)*Ap_mat(2,2)+...
            Ap_mat(1,2)*Ap_mat(2,1)-1);
        ddet = dAp(3,3)*Ap_mat(1,1)+Ap_mat(3,3)*dAp(1,1)+dAp(3,3)*Ap_mat(2,2)+...
            Ap_mat(3,3)*dAp(2,2)-dAp(3,3)*Ap_mat(1,1)*Ap_mat(2,2)-Ap_mat(3,3)*...
            dAp(1,1)*Ap_mat(2,2)-Ap_mat(3,3)*Ap_mat(1,1)*dAp(2,2)+dAp(3,3)*...
            Ap_mat(1,2)*Ap_mat(2,1)+Ap_mat(3,3)*dAp(1,2)*Ap_mat(2,1)+Ap_mat(3,3)*...
            Ap_mat(1,2)*dAp(2,1)-dAp(3,3)-dAp(1,1)-dAp(2,2)+dAp(1,1)*...
            Ap_mat(2,2)+Ap_mat(1,1)*dAp(2,2)-dAp(1,2)*Ap_mat(2,1)-Ap_mat(1,2)*...
            dAp(2,1);
        Am11 = (Ap_mat(3,3)-1)*(Ap_mat(2,2)-1);
        Am12 = -Ap_mat(1,2)*(Ap_mat(3,3)-1);
        Am13 = Ap_mat(1,3)+Ap_mat(1,2)*Ap_mat(2,3)-Ap_mat(1,3)*Ap_mat(2,2);
        Am21 = -Ap_mat(2,1)*(Ap_mat(3,3)-1);
        Am22 = (Ap_mat(1,1)-1)*(Ap_mat(3,3)-1);
        Am23 = Ap_mat(2,3)-Ap_mat(1,1)*Ap_mat(2,3)+Ap_mat(1,3)*Ap_mat(2,1);
        Am33 = -Ap_mat(1,1)-Ap_mat(2,2)+Ap_mat(1,1)*Ap_mat(2,2)-Ap_mat(1,2)*...
            Ap_mat(2,1)+1;
        dAm11 = dAp(3,3)*Ap_mat(2,2)+Ap_mat(3,3)*dAp(2,2)-dAp(3,3)-dAp(2,2);
        dAm12 = -dAp(1,2)*Ap_mat(3,3)-Ap_mat(1,2)*dAp(3,3)+dAp(1,2);
        dAm13 = dAp(1,3)+dAp(1,2)*Ap_mat(2,3)+Ap_mat(1,2)*dAp(2,3)-dAp(1,3)*...
            Ap_mat(2,2)-Ap_mat(1,3)*dAp(2,2);
        dAm21 = -dAp(2,1)*Ap_mat(3,3)-Ap_mat(2,1)*dAp(3,3)+dAp(2,1);
        dAm22 = dAp(1,1)*Ap_mat(3,3)+Ap_mat(1,1)*dAp(3,3)-dAp(1,1)-dAp(3,3);
        dAm23 = dAp(2,3)-dAp(1,1)*Ap_mat(2,3)-Ap_mat(1,1)*dAp(2,3)+dAp(1,3)*...
            Ap_mat(2,1)+Ap_mat(1,3)*dAp(2,1);
        dAm33 = -dAp(1,1)-dAp(2,2)+dAp(1,1)*Ap_mat(2,2)+Ap_mat(1,1)*dAp(2,2)...
            -dAp(1,2)*Ap_mat(2,1)-Ap_mat(1,2)*dAp(2,1);
        Am = [Am11,Am12,Am13;Am21,Am22,Am23;0,0,Am33];
        dAm = [dAm11,dAm12,dAm13;dAm21,dAm22,dAm23;0,0,dAm33];
        
        dInvA = 1/(det)^2*(-1)*ddet*Am+1/det*dAm;
        dFCM_dmu(1:end,1:end,i) = dPhi_dmu(3+3*lenHarm+1+1:end,1:3,i)*...
            (inv(eye(3,3)-Ap_mat))*Np_mat+Hp*dInvA*Np_mat+...
            Hp*inv(eye(3,3)-Ap_mat)*dPhi_dmu(1:3,4:3*lenHarm+4,i)+...
            dPhi_dmu(3*lenHarm+5:end,4:4+3*lenHarm,i);
        Ctemp = dFCM_dmu(1:end,1:end,i)*Zi_arr(1:end,1);
        Ctemp = Ctemp(2*lenHarm+1);
        J(7,i) = (Vdc+2*Rdc*Idc)*Ctemp;
    end
    
    J(7,7) = dIdc_dalpha*(Vdc+2*Rdc*Idc);
    %******************
    %At this point the Matrices should be solved for and I should be able
    %to immediatly obtain M1, M2, M3
    
    M = zeros(7,1);
    for i = 1:6
        At = Mt_arr(1:2,1:2,i);
        Nt = Mt_arr(1:2,4:end,i);
        M(i) = 0-Mi(i,1:end)*(At*inv(ThetaM-Ap(1:2,1:2,i))*Np(1:2,1:end,i)+Nt)*Zi_arr(1:end,i);
    end
    M(7) = Pload - (Idc*Vdc+Rdc*Idc^2);
    
    new_mu = mu + J\M;
    
    flag2 = 0;
    for i = 1:7
        if abs(M(i)) > tolerance
            flag2 = 1;
            break;
        end
    end
    
    if flag2 == 0
        break;
    end

    mu = new_mu;
    alpha = new_mu(7);
    
    if max(mu(1:6)) > pi/2
        flag_error = flag_error +1;
    else
        flag_error = 0;
    end
    if(flag_error > 1)
        error = 1;
        break;
    end
    flag2 = 0;
    countIter = countIter + 1;
end

alpha = mu(7);
%alpha = 0;
RotationMatrix2 = zeros(length(FCM));
count = 1;
for i  = -h:1:h
    RotationMatrix2(count:count+1,count:count+1) = [cos(-i*(alpha-phase)),...
        -sin(-i*(alpha-phase));sin(-i*(alpha-phase)) cos(-i*(alpha-phase))];
    count = count + 2;
end
for i  = 0:1:h
    RotationMatrix2(count:count+1,count:count+1) = [cos(-i*(alpha-phase)),...
        -sin(-i*(alpha-phase));sin(-i*(alpha-phase)) cos(-i*(alpha-phase))];
    count = count + 2;
end
%Not sure if this is the correct rotation matrix...wouldn't the rotation
%matrix be the inverse of the RotationMatrix(1:end,1:end,1) matrix???
tempM = eye(3*lenHarm+1);
tempM(2*lenHarm+3:end,2*lenHarm+1:end) = 2*tempM(2*lenHarm+3:end,2*lenHarm+1:end);
I = tempM*RotationMatrix2*I;
mu_out = mu(1:6);
%FCM = tempM*RotationMatrix2*FCM;

[Ia Ib Ic fundCurrents] = getI_abc_fundCurrents(I,h);
