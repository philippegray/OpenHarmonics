function J12_entry = calcJ12_6(Ap,Np,Mt,dMp,Zi_arr,Zi,Mi,dRotation)
%global lenHarm

At = Mt(1:3,1:3);
Nt = Mt(1:3,4:end);
%dAt = zeros(3,3);
%dNt = zeros(3,3*lenHarm+1);
dAp = dMp(1:3,1:3);
dNp = dMp(1:3,4:end);

dInvImMtx = get_dInvImMtx_6(Ap,dAp);
InvImMtx = inv(eye(3)-Ap);

J12_entry = Mi*(At*InvImMtx*dNp+At*dInvImMtx*Np)*Zi_arr+Mi*(At*...
    InvImMtx*Np+Nt)*dRotation*Zi;