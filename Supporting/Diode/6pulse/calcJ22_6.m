%The inputs are:
%1. Ap(1:end,1:end,i)
%2. Np(1:end,1:end,i),
%3. Mt_arr(1:end,1:end,i,2)
%4. derivMt(1:end,1:end,i,3)
%5. dMp_dmu(1:3,1:3,j,i,1)
%6. Zi_arr(1:end,i)
%7. Mi(i,1:end),i-j);

function J22_entry = calcJ22_6(Ap,Np,dMp,Zi_arr,Zi,Mi,dRotation,constMtx)
global Lout Rout L R

dAp = dMp(1:3,1:3);
dNp = dMp(1:3,4:end);
        
dInvImMtx = get_dInvImMtx_6(Ap,dAp);
InvImMtx = inv(eye(3)-Ap);

J22_entry = (L*(2*R+Rout)/(2*L+Lout)-R)*Mi*((dInvImMtx*Np+InvImMtx*...
    dNp)*Zi_arr+InvImMtx*Np*dRotation*Zi)+constMtx*dRotation*Zi;
