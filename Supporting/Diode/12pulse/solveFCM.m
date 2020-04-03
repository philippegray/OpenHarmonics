function [FCM expmCondMuPi_3 expmCondMu expmCommMu Hp Ap_mat Np_mat Qp]=...
    solveFCM(lenHarm,Ahat,Acomm,Ncomm,Acond,Ncond,mu,gamma)

arrLen = 3+3*lenHarm+1+3*lenHarm+1;

expmCommMu = zeros(arrLen,arrLen,12);
expmCondMu = zeros(arrLen,arrLen,12);
expmCondMuPi_3 = zeros(arrLen,arrLen,12);

Ahat_mod = Ahat;
Phi = eye(arrLen,arrLen);
for i = 1:12
    Ahat_mod(1:3,1:3+3*lenHarm+1) = [Acomm(1:end,1:end,i),Ncomm(1:end,1:end,i)];
    expmCommMu(1:end,1:end,i) = expm(Ahat_mod*mu(i));
    Ahat_mod(1:3,1:3+3*lenHarm+1) = [Acond(1:end,1:end,i),Ncond(1:end,1:end,i)];
    expmCondMu(1:end,1:end,i) = expm(-Ahat_mod*mu(i));
    expmCondMuPi_3(1:end,1:end,i) = expm(Ahat_mod*gamma(i));
    
    Phi = expmCondMuPi_3(1:end,1:end,i)*...
        expmCondMu(1:end,1:end,i)*expmCommMu(1:end,1:end,i)*Phi;
end

Hp = Phi(4+3*lenHarm+1:end,1:3);
Ap_mat = Phi(1:3,1:3);
Np_mat = Phi(1:3,4:3*lenHarm+4);
Qp = Phi(3*lenHarm+4+1:end,4:3*lenHarm+4);

FCM = Hp*(inv(eye(3,3)-Ap_mat))*Np_mat+Qp;
