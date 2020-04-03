function J = calcJ(Mi,Vab_coeff,Mt,dMt_dmu,dMpmu,dMpgamma,Ap,Np,...
    Zi_arr,Zi,Ts,lenHarm,L,R,Ldc,Rdc,derivRotation)

J11 = zeros(12,12);
J12 = zeros(12,13);
J21 = zeros(13,12);
J22 = zeros(13,13);

constMtx = zeros(3,3*lenHarm+1);
for i = 1:2
    for j = i:2:2*lenHarm
        constMtx(i,j) = 1;
    end
end
for i = 1:2:lenHarm+1
    constMtx(3,2*lenHarm+i) = 1;
end

for i = 1:12
    for j = 1:12
        %Entries in Quad. (1,1)
        J11(i,j) = calcJ11(lenHarm,Ap(1:end,1:end,i),Np(1:end,1:end,i),...
            Mt(1:end,1:end,i),dMt_dmu(1:end,1:end,i),...
            dMpmu(1:end,1:end,j,i),Zi_arr(1:end,i),Mi(i,1:end),i-j);
        if j < i
            condMtx = derivRotation(1:end,1:end,i);
        else
            condMtx = zeros(3*lenHarm+1,3*lenHarm+1);
        end
        %Entries in Quad. (1,2)
        J12(i,j) = calcJ12(Ap(1:end,1:end,i),Np(1:end,1:end,i),...
            Mt(1:end,1:end,i),dMpgamma(1:end,1:end,j,i),...
            Zi_arr(1:end,i),Zi,Mi(i,1:end),condMtx);
    end
    condMtx = derivRotation(1:end,1:end,i);
    J12(i,13) = calcJ13(Mt(1:end,1:end,i),Ap(1:end,1:end,i),Np(1:end,1:end,i),...
        Zi,Mi(i,1:end),condMtx);
end
for i = 1:13
    for j = 1:12
        %Entries in Quad. (2,1)
        J21(i,j) = calcJ21(Ap(1:end,1:end,Ts(i)),Np(1:end,1:end,Ts(i)),... %might be looking @wrong value??
            dMpmu(1:end,1:end,j,Ts(i)),Zi_arr(1:end,i),L,R,Ldc,Rdc);
        if j < i
            condMtx = derivRotation(1:end,1:end,i);
        else
            condMtx = zeros(3*lenHarm+1,3*lenHarm+1);
        end
        J22(i,j) = calcJ22(Ap(1:end,1:end,Ts(i)),Np(1:end,1:end,Ts(i)),...
            dMpgamma(1:end,1:end,j,Ts(i)),Zi_arr(1:end,i),Zi,...
            condMtx,[Vab_coeff(Ts(i),1:end),L/(4*L+Ldc)]*constMtx,L,R,Ldc,Rdc);
    end
    condMtx = derivRotation(1:end,1:end,i);
    J22(i,13) = calcJ23(Ap(1:end,1:end,Ts(i)),Np(1:end,1:end,Ts(i)),...
        Zi,condMtx,[Vab_coeff(Ts(i),1:end),L/(4*L+Ldc)]*constMtx,L,R,Ldc,Rdc);
end

J = [J11,J12;J21,J22];
