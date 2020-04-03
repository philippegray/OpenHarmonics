function J23_entry = calcJ23_6(Ap,Np,Zi,Mi,dRotation,constMtx)
global Lout Rout L R

InvImMtx = inv(eye(3)-Ap);

J23_entry = (L*(2*R+Rout)/(2*L+Lout)-R)*Mi*InvImMtx*Np*(-dRotation)*Zi+...
    constMtx*(-dRotation)*Zi;