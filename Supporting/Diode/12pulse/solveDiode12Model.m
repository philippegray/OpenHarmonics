function [Ia Ib Ic fundCurrents phaseAng gamma mu Pload error Rdc] = ...
    solveDiode12Model(Valphabeta,L,R,Rdc,Ldc,h,Vdc,Pdc,mu_i,gamma_i,phase_i)

global SymMtx invSymMtx

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
error = 0;
flag_error = 0;
%model inputs
Vout = Vdc;
[w,Vbase,tol,f,Tac,lenHarm,CTF,CTFsm,...
    invCTF,invCTFsm,Mi,phase,Ts,Vab_coeff,Zi] = inputs_12pulse_diode(Valphabeta,L,R,Ldc,Rdc,h,Vdc);

%adding circuit equations*****
[Acomm,Ncomm,Acond,Ncond,Omegat,Ahat,H] = input_circuiteqns(h,R,L,Rdc,...
    Ldc,lenHarm,CTFsm,invCTFsm,Tac);
%*****************************

mu = mu_i;
phase = phase_i;
gamma = gamma_i;

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
    %     Jnew = zeros(13,13);
    %     Jnew(1:6,1:6) = J(1:6,1:6);
    %     Jnew(7:12,7:12) = J(13:18,13:18);
    %     Jnew(13,13) = J(25,25);
    %     Jnew(1:6,7:12) = J(1:6,13:18);
    %     Jnew(1:6,13) = J(1:6,25);
    %     Jnew(7:12,1:6) = J(13:18,1:6);
    %     Jnew(7:12,13) = J(13:18,25);
    %     Jnew(13,1:6) = J(25,1:6);
    %     Jnew(13,7:12) = J(25,13:18);
    %     Mnew = [M(1:6);M(13:18);M(25)];
    %
    %     Idc = I(2*lenHarm+1);
    %     Vdc = Zi(2*lenHarm+1);
    %     constraintVar = [mu(1:6);gamma(1:6);phase];
    %     constraintVar = constraintVar + Jnew\Mnew;
    %     mu = [constraintVar(1:6);constraintVar(1:6)];
    %     gamma = [constraintVar(7:12);constraintVar(7:12)];
    %     phase = constraintVar(13);
    
    
    Idc = I(2*lenHarm+1);
    if count == 0
        Iinit = Idc;
    end
    
    Vdc = Zi(2*lenHarm+1);
    constraintVar = [mu;gamma;phase];
    constraintVar = constraintVar + inv(J)*M;
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
    if(flag_error >= 100)
        break;
    else
        if max(mu(1:12)) > pi/3 || min(mu(1:12)) < -pi/3
            flag_error = flag_error +1;
            if max(mu(1:6)) > pi || min(mu(1:6)) < -pi
                flag_error = flag_error+3;
            end
        else
            flag_error = 0;
        end
        if(flag_error > 2)
            error = 1;
            gamma = gamma_i;
            phase = phase_i;
            mu = mu_i;
            flag_error = 100;
            break;
        end
    end
    
    %
    %     for i = 1:12
    %         if mu(i) > pi || mu(i) < 0
    %             mu(i) = pi/6;
    %         elseif mu(i) < 0
    %             mu(i) = 0;
    %         end
    %     end
    
    count = count + 1;
    if count > 25
        flag_error = 100;
    end
end

for i = 1:12
   if mu(i) < 0 || mu(i) > gamma(i)
       error = 1;
   end
end
if error ~= 0
    Idc = Iinit;
end

phaseAng = phase;
Pload = Rdc*Idc^2;
[Ia Ib Ic fundCurrents] = getI_abc_fundCurrents(I,h);
