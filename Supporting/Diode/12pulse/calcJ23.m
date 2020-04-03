function J23_entry = calcJ23(Ap,Np,Zi,dRotation,constMtx,L,R,Ldc,Rdc)

InvImMtx = inv(eye(3)-Ap);

J23_entry = (L*(4*R+Rdc)/(4*L+Ldc)-R)*[0 0 1]*InvImMtx*Np*(-dRotation)*Zi+...
    constMtx*(-dRotation)*Zi;