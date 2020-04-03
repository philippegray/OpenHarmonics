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

tic;
%model inputs
[w,h,R,L,Rdc,Ldc,Vbase,tol,f,Tac,Vout,lenHarm,CTF,CTFsm,...
    invCTF,invCTFsm,Mi,phase,Ts,Zi,Vab_coeff] = inputs_12pulse_diode();

%adding circuit equations*****
[Acomm,Ncomm,Acond,Ncond,Omegat,Ahat,H] = input_circuiteqns(h,R,L,Rdc,...
    Ldc,lenHarm,CTFsm,invCTFsm,Tac);
%*****************************

mu = (pi/12)*ones(12,1);%(0.432)*ones(13,1); %Initialization
gamma = pi/6*ones(12,1);
count = 0;
while(true)
    
    %solving for the rotated initial state vector
    Zi_arr = rotMtx(lenHarm,h,gamma,phase,Zi); %Fine.
    
    %solving for FCM and exponential matrices
    [FCM expmCondMuPi_3 expmCondMu expmCommMu Hp Ap_mat Np_mat Qp] = ...
        solveFCM(lenHarm,Ahat,Acomm,Ncomm,Acond,Ncond,mu,gamma); %exmpCommMu matrices equate
    
    I = FCM*Zi_arr(1:end,1); %this matches with the thyristor 12 pulse! so FCM is 100% correct.
        
    %Solve Mp
    [Mp dMpmu dMpgamma Ap Np] = solveMp(lenHarm,Acond,Ncond,...
    Acomm,Ncomm,Omegat,expmCondMuPi_3,expmCondMu,expmCommMu,Ts);
    
    %solving Mt and dMt_d(mu)
    [Mt dMt_dmu] = solveMt(lenHarm,Acomm,Ncomm,Omegat,expmCommMu); %Fine
    
    %solving for the derivative of the rotational matrix
    derivRotation = deriv_rotMtx(lenHarm,h,gamma,phase);
    
    %calculate M - not sure if this is fully correct - especially the Mi
    %values.
    M = calcM(Mi,Vab_coeff,Mt,Ap,Np,Zi_arr,Ts,lenHarm,L,R,Ldc,Rdc);
    
    %Works up to here!
    %calculate Jacobian J
    J = calcJ(Mi,Vab_coeff,Mt,dMt_dmu,dMpmu,dMpgamma,Ap,Np,Zi_arr,Zi,...
        Ts,lenHarm,L,R,Ldc,Rdc,derivRotation);
    
    Idc = I(2*lenHarm+1);
    Vdc = Zi(2*lenHarm+1);
    constraintVar = [mu;gamma;phase];
    constraintVar = constraintVar + J\M;
    mu = constraintVar(1:12);
    gamma = constraintVar(13:24);
    phase = constraintVar(25);
    
    flag2 = 0;
    for i = 1:25
        if abs(M(i)) > tol
            flag2 = 1;
            break;
        end
    end
    
    if flag2 == 0
        break;
    end
%     for i = 1:12
%         if mu(i) > pi || mu(i) < 0
%             mu(i) = pi/6;
%         elseif mu(i) < 0
%             mu(i) = 0;
%         end
%     end
    
    if count > 14
        count
        break;
    end

    count = count + 1;
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
