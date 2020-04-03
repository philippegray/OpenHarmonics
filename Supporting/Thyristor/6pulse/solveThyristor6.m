%can apply this exact template for the other thyristor and diode models
function [Ia,Ib,Ic,fundCurrents,alpha,mu] = solveThyristor6(Valphabeta,L,R,...
    Rdc,Ldc,h,Vdc,P_in)

global OpenDSSFileLoc

P = P_in;
errorFile = fopen([OpenDSSFileLoc,'errorFile.txt'], 'w');

%Pmax_est = 1/(Rdc+2*R)*(-Vdc_in^2+3*sqrt(3)/pi*Vdc_in*cos(0));
%Pmin_est = 1/(Rdc+2*R)*(-Vdc_in^2+3*sqrt(3)/pi*Vdc_in*cos(60*pi/180));

if Vdc >= 1.7321 || Vdc < 0
    if Vdc >= 1.7321
        fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Thyristor 6 pulse model will not conduct any current because Vdc_pu > 1.73\n',c(4),c(5),c(6));
    else
        fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Thyristor 6 pulse model will not conduct any current Vdc_pu must be > 0\n',c(4),c(5),c(6));
    end
    Iphase = zeros(h,3);
    
    for i = 1:h
        Iphase(i,1) = i;
    end
    Iphase(1,2) = 100;
    Ia = Iphase;
    Ib = Iphase;
    Ic = Iphase;
    fundCurrents = [0,0;0,-120;0,+120];
    alpha = 0;
    mu = zeros(6,1);
    return;
end

if P > 3 %restricting to 2 pu
    c = clock;
    fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs P > Pmax = 3.0 pu. Setting P = 3.0 pu.\n',c(4),c(5),c(6));
    Iphase = zeros(h,3);
    
    for i = 1:h
        Iphase(i,1) = i;
    end
    Iphase(1,2) = 100;
    Ia = Iphase;
    Ib = Iphase;
    Ic = Iphase;
    fundCurrents = [0,0;0,-120;0,+120];
    alpha = 0;
    mu = zeros(6,1);
    return;
end
if P <= 0 %restricting to 2 pu
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
    alpha = 0;
    mu = zeros(6,1);
    fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs P < Pmin = 0.0 pu.\n',c(4),c(5),c(6));
    return;
end

%Initializations
alpha_i = 0;
alpha = alpha_i;
mu = -100;
Vs_mag = abs(Valphabeta(2*h+3)+1i*Valphabeta(2*h+4));

while (alpha > pi*3/2 || alpha < 0 || mu < 0 || mu > pi/3)
    convergance = 0;
    alpha_i = alpha_i+pi/3/100;
    alpha = alpha_i;
    if alpha > pi*3/2
        c = clock;
        fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs 12 pulse thyristor model cannot find solution. Try another Vdc value.\n',c(4),c(5),c(6));
        Iphase = zeros(h,3);
        
        for i = 1:h
            Iphase(i,1) = i;
        end
        Iphase(1,2) = 100;
        Ia = Iphase;
        Ib = Iphase;
        Ic = Iphase;
        fundCurrents = [0,0;0,-120;0,+120];
        alpha = 0;
        mu = zeros(12,1);
        return;
    end
    while(convergance == 0)
        checkVal = 2*pi/3/sqrt(3)*Vdc/Vs_mag - cos(alpha);
        while(checkVal < -1 || checkVal > 1)
            alpha_i = alpha_i+pi/3/100;
            alpha = alpha_i;
            checkVal = 2*pi/3/sqrt(3)*Vdc/Vs_mag - cos(alpha);
            if alpha > pi*3/2
                c = clock;
                fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs 12 pulse thyristor model cannot find solution. Try another Vdc value.\n',c(4),c(5),c(6));
                Iphase = zeros(h,3);
                
                for i = 1:h
                    Iphase(i,1) = i;
                end
                Iphase(1,2) = 100;
                Ia = Iphase;
                Ib = Iphase;
                Ic = Iphase;
                fundCurrents = [0,0;0,-120;0,+120];
                alpha = 0;
                mu = zeros(12,1);
                return;
            end
        end
        mu = acos(checkVal)-alpha;
        Idc_calc = 1/2/L*(-sqrt(3)*cos(alpha+mu)+sqrt(3)*cos(alpha))*Vs_mag;
        
        J = sqrt(3)/2/L*(-sqrt(3)*sin(alpha)-sin(alpha)/sqrt(1-(pi/3/sqrt(3)*Vdc/Vs_mag-cos(alpha))))*Vs_mag;
        
        M = P-Vdc*Idc_calc;
        
        
        if abs(M) < 0.00001
            convergance = 1;
            break;
        end
        
        alpha = alpha+1/J*M;
    end
end

mu_prev = mu*ones(6,1);
alpha_prev = alpha;

center = 2*h+1;
convergance = 0;
%****************

for i = 1:2:h
    V = Valphabeta(center-i*2:center+i*2+1);
    while(convergance == 0)      
        [Ia Ib Ic fundCurrents alpha mu error] =...
            solveThyristor6Model(V,L,R,Rdc,Ldc,i,0,P,mu_prev,alpha_prev);

        convergance = 1;
        
        mu_prev = mu;
        alpha_prev = alpha;
        
    end
    %mu_prev = mu;
    %alpha_prev = alpha;
    convergance = 0;
end


















% 
% 
% while(error ~= 0)
%     [Ia,Ib,Ic,fundCurrents,alpha,mu_out,error] = ...
%         solveThyristor6Model(Valphabeta,L,R,Rdc,Ldc,h,Vdc,P);
%     
%     if error == 1
%         c = clock;
%         coeff = 0.2;
%         Rdc = Rdc*coeff;
%         fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Changing Rdc from %1.2f to %1.2f\n',c(4),c(5),c(6),Rdc/coeff,Rdc);
%     else
%         for i = 1:6
%             if mu_out(i) < 0
%                 c = clock;
%                 coeff = 2;
%                 Rdc = Rdc*coeff;
%                 fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Changing Rdc from %1.2f to %1.2f\n',c(4),c(5),c(6),Rdc/coeff,Rdc);
%                 error = 2;
%                 break;
%             end
%         end
%     end
% end


fclose(errorFile);