function J13_entry = calcJ13_6(Mt,Ap,Np,Zi,Mi,dRotation)

At = Mt(1:3,1:3);
Nt = Mt(1:3,4:end);
InvImMtx = inv(eye(3)-Ap);

J13_entry = Mi*(At*InvImMtx*Np+Nt)*(-dRotation)*Zi;