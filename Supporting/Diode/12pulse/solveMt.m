function [Mt dMt_dmu] = solveMt(lenHarm,Acomm,Ncomm,Omegat,expmCommMu)

maxPos = 3+3*lenHarm+1;
dMt_dmu = zeros(maxPos,maxPos,12);
Mt = zeros(maxPos,maxPos,12);
Mcomm = zeros(maxPos,maxPos);

for i = 1:12
    %Solving Mt and dMt
    Mcomm(1:3,1:3) = Acomm(1:end,1:end,i);
    Mcomm(1:3,4:end) = Ncomm(1:end,1:end,i);
    Mcomm(4:end,4:end) = Omegat;
    
    Mt(1:end,1:end,i) = expmCommMu(1:maxPos,1:maxPos,i);
    dMt_dmu(1:end,1:end,i) = Mcomm*Mt(1:end,1:end,i);
    %******************
end