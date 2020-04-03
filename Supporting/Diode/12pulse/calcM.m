function M = calcM(Mi,Vab_coeff,Mt,Ap,Np,Zi_arr,Ts,lenHarm,L,R,Ldc,Rdc)

M = zeros(25,1);
for i = 1:12
    At = Mt(1:3,1:3,i);
    Nt = Mt(1:3,4:end,i);
    M(i) = 0-Mi(i,1:end)*(At*inv(eye(3,3)-Ap(1:3,1:3,i))*...
        Np(1:3,1:end,i)+Nt)*Zi_arr(1:end,i);
end

constMtx = zeros(3,3*lenHarm+1);
for i = 1:2
    for j = i:2:2*lenHarm
        constMtx(i,j) = 1;
    end
end
for i = 1:2:lenHarm+1
    constMtx(3,2*lenHarm+i) = 1;
end

for i = 1:13 %need to use Ts b/c I am solving for 7 intervals > num of intervals/period
    Idc = [0 0 1]*(inv(eye(3,3)-Ap(1:3,1:3,Ts(i)))*...
        Np(1:3,1:end,Ts(i)))*Zi_arr(1:end,i);
    
    M(i+12) = 0 - ((L*(4*R+Rdc)/(4*L+Ldc)-R)*Idc+...
        [Vab_coeff(Ts(i),1:end),L/(4*L+Ldc)]*constMtx*Zi_arr(1:end,i));
end %why does Zi_arr have to be Ts(i+1) - shouldn't it be just i..
%calculating M(13)*********
