function [Mp dMpmu dMpgamma Ap Np] = solveMp(lenHarm,Acond,Ncond,...
    Acomm,Ncomm,Omegat,expmCondMuPi_3,expmCondMu,expmCommMu,Ts)

maxPos = 3+3*lenHarm+1;
dMp_dmu = zeros(maxPos,maxPos,12,12); %3rd column contains the partial
dMp_dgamma = zeros(maxPos,maxPos,12,12); %3rd column contains the partial
Mp = zeros(maxPos,maxPos,12);
Ap = zeros(3,3,12);
Np = zeros(3,3*lenHarm+1,12);
dMpgamma = zeros(length(Mp),length(Mp),12,12);
dMpmu = zeros(length(Mp),length(Mp),12,12);
Mcond = zeros(maxPos,maxPos);
Mcomm = zeros(maxPos,maxPos);

for j = 1:12
    Mp(1:end,1:end,j) = eye(maxPos,maxPos);
    for k = 1:12
        dMp_dmu(1:end,1:end,j,k) = eye(maxPos,maxPos);
        dMp_dgamma(1:end,1:end,j,k) = eye(maxPos,maxPos);%zeros(length(derivMp));%
    end
end

%solving for Mp and dMp and Ap and Np
for i = 1:12
    for j = 0:11
        Mcond(1:3,1:3) = Acond(1:end,1:end,Ts(i+j));
        Mcond(1:3,4:end) = Ncond(1:end,1:end,Ts(i+j));
        Mcond(4:end,4:end) = Omegat;
        Mcomm(1:3,1:3) = Acomm(1:end,1:end,Ts(i+j));
        Mcomm(1:3,4:end) = Ncomm(1:end,1:end,Ts(i+j));
        Mcomm(4:end,4:end) = Omegat;
        
        Mp(1:end,1:end,i) = expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
            expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
            expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*Mp(1:end,1:end,i);
        
        for k = 1:12
            if k == (j+1)
                dMp_dmu(1:end,1:end,k,i) = ...
                    expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                    (-Mcond*expmCondMu(1:maxPos,1:maxPos,Ts(i+j))+...
                    expmCondMu(1:maxPos,1:maxPos,Ts(i+j))*Mcomm)*...
                    expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                    dMp_dmu(1:end,1:end,k,i);
                
                dMp_dgamma(1:end,1:end,k,i) = Mcond*...
                    expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                    expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                    expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                    dMp_dgamma(1:end,1:end,k,i);
            else
                dMp_dmu(1:end,1:end,k,i) = ...
                    expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                    expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                    expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                    dMp_dmu(1:end,1:end,k,i);
                
                dMp_dgamma(1:end,1:end,k,i) = ...
                    expmCondMuPi_3(1:maxPos,1:maxPos,Ts(j+i))*...
                    expmCondMu(1:maxPos,1:maxPos,Ts(j+i))*...
                    expmCommMu(1:maxPos,1:maxPos,Ts(i+j))*...
                    dMp_dgamma(1:end,1:end,k,i);
            end
        end
    end
end

%the 4th dimension is whether it be M1,M2,M3 (each has an offset to the initial vector)
%the 3rd dimension is for the derivative with respect to mu1,mu2,mu3

for i = 0:11
    for j = 1:12
        dMpgamma(1:end,1:end,Ts(j+i),i+1) = dMp_dgamma(1:end,1:end,j,i+1);
        dMpmu(1:end,1:end,Ts(j+i),i+1) = dMp_dmu(1:end,1:end,j,i+1);
    end
end
for i = 1:12
    Ap(1:end,1:end,i) = Mp(1:3,1:3,i);
    Np(1:end,1:end,i) = Mp(1:3,4:end,i);
end
