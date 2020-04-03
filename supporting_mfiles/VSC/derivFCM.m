function dFCM = derivFCM(dHp,dNp,dQp,dAp,Hp,Np,Ap)

dInvImMtx = get_dInvImMtx(Ap,dAp);
InvImMtx = inv(eye(3)-Ap);

dFCM = dHp*InvImMtx*Np+Hp*dInvImMtx*Np+Hp*InvImMtx*dNp+dQp;