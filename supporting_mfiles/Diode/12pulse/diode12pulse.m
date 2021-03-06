%Model for Thyristor 12 Pulse Rectifier*****************
%Model takes inputs and produces outputs
%Inputs: U matrix = [Vsab^-h ... Vsab^+h Vdc^0 ... Vdc^+h]
%Outputs: Y vector = [Isab^-h ... Isab^+h Idc^0 ... Idc^+h]
%There are one iterative loops in the program:
%outer loop finds the mu1...mu12 alpha values that gives you the desired
%real power of the converter - one control variable (alpha) = one thing
%controlled (real power)
%in that loop we find an FCM and a jacobian for the outer control loop
%Philippe A. Gray
%April 26, 2012
%********************************************************
cd('C:\Users\Phil\Dropbox\Masters\Matlab\Diode12Pulse\');
clear all;

tic;
%model inputs
[w,h,R,L,Rdc,Ldc,Vbase,tol,f,Tac,Vout,lenHarm,CTF,CTFsm,...
    invCTF,invCTFsm,Mi,phase,Ts,Zi,Vab_coeff] = inputs_12pulse_diode();

%adding circuit equations
[Acomm,Ncomm,Acond,Ncond,Omegat,Ahat,H] = input_circuiteqns(h,R,L,Rdc,...
    Ldc,lenHarm,CTFsm,invCTFsm,Tac);
%************************

mu = (pi/6)*ones(12,1);%(0.432)*ones(13,1); %Initialization
gamma = pi/6*ones(12,1);

while(true)
    
    %solving for the rotated initial state vector
    Zi_arr = rotMtx(lenHarm,h,gamma,phase,Zi); %Fine.
    
    %solving for FCM and exponential matrices
    [FCM expmCondMuPi_3 expmCondMu expmCommMu Hp Ap_mat Np_mat Qp] = ...
        solveFCM(lenHarm,Ahat,Acomm,Ncomm,Acond,Ncond,mu,gamma); %exmpCommMu matrices equate
    
    I = FCM*Zi_arr(1:end,1); %this matches with the thyristor 12 pulse! so FCM is 100% correct.
        
    %test
    
    maxPos = 3+3*lenHarm+1;
    dMp_dmu = zeros(maxPos,maxPos,12,12); %3rd column contains the partial
    dMp_dgamma = zeros(maxPos,maxPos,12,12); %3rd column contains the partial
    Mp = zeros(maxPos,maxPos,12);
    Ap = zeros(3,3,12);
    Np = zeros(3,3*lenHarm+1,12);
    dMpgamma = zeros(length(Mp),length(Mp),12,12);
    dMpmu = zeros(length(Mp),length(Mp),12,12);
    Mcond = zeros(maxPos,maxPos);
    Mcomm = zeros(maxPos,maxPos);
    
    for j = 1:12
        Mp(1:end,1:end,j) = eye(maxPos,maxPos);
        for k = 1:12
            dMp_dmu(1:end,1:end,j,k) = eye(maxPos,maxPos);
            dMp_dgamma(1:end,1:end,j,k) = eye(maxPos,maxPos);%zeros(length(derivMp));%
        end
    end
    
    for i = 1:12
        for j = 0:11
            Mcond(1:3,1:3) = Acond(1:end,1:end,Ts(i+j));
            Mcond(1:3,4:end) = Ncond(1:end,1:end,Ts(i+j));
            Mcond(4:end,4:end) = Omegat;
            Mcomm(1:3,1:3) = Acomm(1:end,1:end,Ts(i+j));
            Mcomm(1:3,4:end) = Ncomm(1:end,1:end,Ts(i+j));
            Mcomm(4:end,4:end) = Omegat;
            
            Mp(1:end,1:end,i) = expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*Mp(1:end,1:end,i);
            
            for k = 1:12
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
    
    for i = 0:11
        for j = 1:12
            dMpgamma(1:end,1:end,Ts(j+i),i+1) = dMp_dgamma(1:end,1:end,j,i+1);
            dMpmu(1:end,1:end,Ts(j+i),i+1) = dMp_dmu(1:end,1:end,j,i+1);
        end
    end
    for i = 1:12
        Ap(1:end,1:end,i) = Mp(1:3,1:3,i);
        Np(1:end,1:end,i) = Mp(1:3,4:end,i);
    end

    %***********
    
    %solving Mp and dMp_d(mu,gamma)
    %[Mp dMp_dmu dMp_dgamma Ap Np] = solveMp(lenHarm,Acond,Ncond,Acomm,...
    %    Ncomm,Omegat,expmCondMuPi_3,expmCondMu,expmCommMu,Ts); %the exponential matrices are all correct...they match with equivalent thyristor model.
    
    %solving Mt and dMt_d(mu)
    [Mt dMt_dmu] = solveMt(lenHarm,Acomm,Ncomm,Omegat,expmCommMu); %Fine
    
    %solving for the derivative of the rotational matrix
    derivRotation = deriv_rotMtx(lenHarm,h,gamma,phase);
    
    %calculate M - not sure if this is fully correct - especially the Mi
    %values.
    M = calcM(Mi,Vab_coeff,Mt,Ap,Np,Zi_arr,Ts,lenHarm,L,R,Ldc,Rdc);
    
    %calculate Jacobian J
    
    
    %**************************************************
    J = zeros(13,13);
    derivInv2 = zeros(3,3);
    for i = 1:12
        for j = 1:12
            At = Mt(1:3,1:3,i);
            if i == j
                dAt = dMt(1:3,1:3,i);
                dNt = dMt(1:3,4:end,i);
            else
                dAt = zeros(3,3);
                dNt = zeros(3,3*lenHarm+1);
            end
            dAp = dMp_dmu(1:3,1:3,j,i);
            dNp = dMp_dmu(1:3,4:end,j,i);
            invA = inv(ThetaM-Ap(1:3,1:3,i));
            %FIX THIS!!!!!!**********************
            det = (1-Ap(3,3,i))*(1-Ap(1,1,i)-Ap(2,2,i)+Ap(1,1,i)*...
                Ap(2,2,i)-Ap(1,2,i)*Ap(2,1,i));

            ddet_mu = (-1+Ap(2,2,i)+Ap(3,3,i)-Ap(3,3,i)*Ap(2,2,i))*dAp(1,1)+...
                (-Ap(2,1,i)+Ap(3,3,i)*Ap(2,1,i))*dAp(1,2)+...
                (-Ap(1,2,i)+Ap(3,3,i)*Ap(1,2,i))*dAp(2,1)+...
                (-1+Ap(1,1,i)+Ap(3,3,i)-Ap(3,3,i)*Ap(1,1,i))*dAp(2,2)+...
                (-1+Ap(1,1,i)+Ap(2,2,i)-Ap(1,1,i)*Ap(2,2,i)+Ap(1,2,i)*Ap(2,1,i))*dAp(3,3);
            derivInv2(1,1) = det*(-dAp(2,2)-dAp(3,3)+Ap(2,2,i)*dAp(3,3)+...
                Ap(3,3,i)*dAp(2,2))-ddet_mu*(1-Ap(2,2,i)-Ap(3,3,i)+Ap(2,2,i)*Ap(3,3,i));
            derivInv2(1,2) = det*(dAp(1,2)-Ap(1,2,i)*dAp(3,3)-Ap(3,3,i)*dAp(1,2))-...
                ddet_mu*(Ap(1,2,i)-Ap(1,2,i)*Ap(3,3,i));
            derivInv2(2,1) = det*(dAp(2,1)-Ap(2,1,i)*dAp(3,3)-Ap(3,3,i)*dAp(2,1))-...
                ddet_mu*(Ap(2,1,i)-Ap(2,1,i)*Ap(3,3,i));
            derivInv2(2,2) = det*(-dAp(1,1)-dAp(3,3)+Ap(1,1,i)*dAp(3,3)+...
                Ap(3,3,i)*dAp(1,1))-ddet_mu*(1-Ap(1,1,i)-Ap(3,3,i)+Ap(1,1,i)*Ap(3,3,i));
            derivInv2(3,3) = det*(-dAp(1,1)-dAp(2,2)+Ap(1,1,i)*dAp(2,2)+...
                Ap(2,2,i)*dAp(1,1)-Ap(1,2,i)*dAp(2,1)-Ap(2,1,i)*dAp(1,2))-...
                ddet_mu*(1-Ap(1,1,i)-Ap(2,2,i)+Ap(1,1,i)*Ap(2,2,i)-Ap(1,2,i)*Ap(2,1,i));
            derivInv2 = 1/(det)^2*derivInv2;
            %derivInv = -ddet_dmu/det^2*AM +1/det*dAM_dmu;
            J(i,j) = Mi_deriv(i,1:end)*(dAt*invA*Np(1:3,1:end,i)+At*invA*dNp+...
                At*derivInv2*Np(1:3,1:end,i)+dNt)*Zi_arr(1:end,i);
            %***************************************
        end
    end
      
    for i = 1:12
        At = Mt(1:3,1:3,i);
        Nt = Mt(1:3,4:end,i);
        J(i,13) = Mi_deriv(i,1:end)*(At*inv(ThetaM-Ap(1:3,1:3,i))*...
            Np(1:3,1:end,i)+Nt)*derivRotation(1:end,1:end,i)*Zi;
    end
    
    dFCM_dmu = zeros(length(FCM),length(FCM),12);
    derivInv2 = zeros(3,3);
    derivAp = zeros(3,3);
    det = (1-Ap_mat(3,3))*(1-Ap_mat(1,1)-Ap_mat(2,2)+Ap_mat(1,1)*Ap_mat(2,2)-...
        Ap_mat(1,2)*Ap_mat(2,1));
    for i = 1:12
        derivAp = dPhi_dmu(1:3,1:3,i);
        ddet_mu = (-1+Ap_mat(2,2)+Ap_mat(3,3)-Ap_mat(3,3)*Ap_mat(2,2))*derivAp(1,1)+...
            (-Ap_mat(2,1)+Ap_mat(3,3)*Ap_mat(2,1))*derivAp(1,2)+...
            (-Ap_mat(1,2)+Ap_mat(3,3)*Ap_mat(1,2))*derivAp(2,1)+...
            (-1+Ap_mat(1,1)+Ap_mat(3,3)-Ap_mat(3,3)*Ap_mat(1,1))*derivAp(2,2)+...
            (-1+Ap_mat(1,1)+Ap_mat(2,2)-Ap_mat(1,1)*Ap_mat(2,2)+Ap_mat(1,2)*Ap_mat(2,1))*derivAp(3,3);
        derivInv2(1,1) = det*(-derivAp(2,2)-derivAp(3,3)+Ap_mat(2,2)*derivAp(3,3)+...
            Ap_mat(3,3)*derivAp(2,2))-ddet_mu*(1-Ap_mat(2,2)-Ap_mat(3,3)+Ap_mat(2,2)*Ap_mat(3,3));
        derivInv2(1,2) = det*(derivAp(1,2)-Ap_mat(1,2)*derivAp(3,3)-Ap_mat(3,3)*derivAp(1,2))-...
            ddet_mu*(Ap_mat(1,2)-Ap_mat(1,2)*Ap_mat(3,3));
        derivInv2(2,1) = det*(derivAp(2,1)-Ap_mat(2,1)*derivAp(3,3)-Ap_mat(3,3)*derivAp(2,1))-...
            ddet_mu*(Ap_mat(2,1)-Ap_mat(2,1)*Ap_mat(3,3));
        derivInv2(2,2) = det*(-derivAp(1,1)-derivAp(3,3)+Ap_mat(1,1)*derivAp(3,3)+...
            Ap_mat(3,3)*derivAp(1,1))-ddet_mu*(1-Ap_mat(1,1)-Ap_mat(3,3)+Ap_mat(1,1)*Ap_mat(3,3));  
        derivInv2(3,3) = det*(-derivAp(1,1)-derivAp(2,2)+Ap_mat(1,1)*derivAp(2,2)+...
            Ap_mat(2,2)*derivAp(1,1)-Ap_mat(1,2)*derivAp(2,1)-Ap_mat(2,1)*derivAp(1,2))-...
            ddet_mu*(1-Ap_mat(1,1)-Ap_mat(2,2)+Ap_mat(1,1)*Ap_mat(2,2)-Ap_mat(1,2)*Ap_mat(2,1));
        derivInv2 = 1/(det)^2*derivInv2;
        dFCM_dmu(1:end,1:end,i) = dPhi_dmu(3+3*lenHarm+1+1:end,1:3,i)*...
            (inv(eye(3,3)-Ap_mat))*Np_mat+Hp*derivInv2*Np_mat+...
            Hp*inv(eye(3,3)-Ap_mat)*dPhi_dmu(1:3,4:3*lenHarm+4,i)+...
            dPhi_dmu(3*lenHarm+5:end,4:4+3*lenHarm,i);
        Ctemp = dFCM_dmu(1:end,1:end,i)*Zi_arr(1:end,1);
        Ctemp = Ctemp(2*lenHarm+1);
        J(13,i) = Vdc*Ctemp;
    end
    
    J(13,13) = dIdc_dalpha*Vdc;
    %******************
    %At this point the Matrices should be solved for and I should be able
    %to immediatly obtain M1, M2, M3
    

    Idc = I(2*lenHarm+1);
    Vdc = Zi(2*lenHarm+1);
    
    %M(7) = (new_mu(7)-mu(7))/mu(7)*100;
    
    new_mu = mu + J\M;
    
    flag2 = 0;
    for i = 1:12
        if abs(M(i)) > tolerance
            flag2 = 1;
            break;
        end
    end
    
    if flag2 == 0
        break;
    end
    for i = 1:12
        if new_mu(i) > pi/6 || new_mu(i) < 0
            new_mu(i) = pi/6;
        elseif new_mu(i) < 0
            new_mu(i) = 0;
        end
    end
    
    if countIter > 19
        countIter
        break;
    end
    %new_mu(7) = 0.0;
    %alpha = 0.3;
    %alpha = 0;
    mu = new_mu;
    alpha = new_mu(13);
    flag2 = 0;
    countIter = countIter + 1;
end
%alpha = mu(7);
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

m = h;
Imag = zeros(lenHarm+1+m,2);
for i = 0:2*h+1-1
   Imag(i+1,1:2) = [i-h,abs(I(2*i+1)+1i*I(2*i+2))];
end
count = 1;
for i = 1:m+1
   Imag(2*h+1+i,1:2) = [i-1,abs(I(2*lenHarm+count)+1i*I(2*lenHarm+count+1))];
   count = count + 2;
end

toc;
