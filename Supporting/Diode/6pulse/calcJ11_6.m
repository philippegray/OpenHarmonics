%The inputs are:
%1. Ap(1:end,1:end,i)
%2. Np(1:end,1:end,i),
%3. Mt_arr(1:end,1:end,i,2)
%4. derivMt(1:end,1:end,i,3)
%5. dMp_dmu(1:3,1:3,j,i,1)
%6. Zi_arr(1:end,i)
%7. Mi(i,1:end),i-j);

function J11_entry = calcJ11_6(Ap,Np,Mt,dMt,dMp,Zi_arr,Mi,imj)
global lenHarm

At = Mt(1:3,1:3);
if imj == 0
    dAt = dMt(1:3,1:3);
    dNt = dMt(1:3,4:end);
else
    dAt = zeros(3,3);
    dNt = zeros(3,3*lenHarm+1);
end
dAp = dMp(1:3,1:3);
dNp = dMp(1:3,4:end);
        
dInvImMtx = get_dInvImMtx_6(Ap,dAp);
InvImMtx = inv(eye(3)-Ap);

J11_entry = Mi*(dAt*InvImMtx*Np+At*InvImMtx*dNp+At*dInvImMtx*Np+dNt)*Zi_arr;

