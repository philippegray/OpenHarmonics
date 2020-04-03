function [Ia Ib Ic fundCurrents ma sigma Rdc error] = solveVSC(Valphabeta,...
                    L,R,mf,Rdc_in,Cdc,Idc_in,h,Vdc_in,Q)

global OpenDSSFileLoc

error = 0;
Vdc = Vdc_in;
Rdc = Rdc_in;
Idc = Idc_in;
% if Idc ~= 0
%     P = Vdc_in/Idc;
% else
%     P = 0;
% end

errorFile = fopen([OpenDSSFileLoc,'errorFile.txt'], 'w');
if Idc > 3 %restricting to 2 pu
   c = clock;
   Idc = 3;
   fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs P > Pmax = 3.0 pu. Setting P = 3.0 pu.\n',c(4),c(5),c(6));
elseif Idc < -3 %restricting to 2 pu
   c = clock;
   Idc = -3;
   fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs P < Pmin = -3.0 pu. Setting P = -3.0 pu.\n',c(4),c(5),c(6));
end
if Vdc <= 0 %restricting to 2 pu
    c = clock;
    Iphase = zeros(h,3);
    for i = 1:h
        Iphase(i,1) = i;
    end
    Iphase(1,2) = 100;
    Ia = Iphase;
    Ib = Iphase;
    Ic = Iphase;
    fundCurrents = [0,0;0,-120;0,+120];
    sigma = 0;
    ma = 0;
    error = 0;
    fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs The Vdc for VSC must be > 0\n',c(4),c(5),c(6));
    return;
end
if Vdc >= 3 %restricting to 2 pu
    c = clock;
    Vdc = 3;
    fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs The Vdc for VSC must be < 3 pu. Setting Vdc to 3.0 pu.\n',c(4),c(5),c(6));
end
if Idc == 0 %restricting to 2 pu
    Idc = 0.001;
end
%     Iphase = zeros(h,3);
%     for i = 1:h
%         Iphase(i,1) = i;
%     end
%     Iphase(1,2) = 100;
%     Ia = Iphase;
%     Ib = Iphase;
%     Ic = Iphase;
%     fundCurrents = [0,0;0,-120;0,+120];
%     sigma = 0;
%     ma = 0;
%     return;
% end
[Ia Ib Ic fundCurrents ma sigma Rdc error] = VSCmodel(Valphabeta,L,R,mf,Rdc,Cdc,Idc,h,Vdc,Q);

if error == 1
    c = clock;
    Iphase = zeros(h,3);
    for i = 1:h
        Iphase(i,1) = i;
    end
    Iphase(1,2) = 100;
    Ia = Iphase;
    Ib = Iphase;
    Ic = Iphase;
    fundCurrents = [0,0;0,-120;0,+120];
    sigma = 0;
    ma = 0;
    return;
end
fclose(errorFile);