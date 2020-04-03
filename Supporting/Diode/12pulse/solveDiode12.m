function [Ia Ib Ic fundCurrents phaseAng gamma mu Pload Rdc] = ...
    solveDiode12(Valphabeta,L,R,Rdc_in,Ldc,h,Vdc_in,P_in,mu_in,gamma_in,phaseAng_in)

global OpenDSSFileLoc

error = 1;
Vdc = Vdc_in;
Rdc = Rdc_in;
P = P_in;
errorFile = fopen([OpenDSSFileLoc,'errorFile.txt'], 'w');
if P_in > 3 %restricting to 2 pu
   c = clock;
   P = 3;
   fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs P > Pmax = 3.0 pu. Setting P = 3.0 pu.\n',c(4),c(5),c(6));
end
if P_in <= 0 %restricting to 2 pu
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
    phaseAng = 0;
    gamma = 0;
    mu = 0;
    Pload = 0;
    return;
end

%Initializations

if mu_in == -1
    Idc_ref = (-Vdc+sqrt(Vdc+4*Rdc*P))/2/Rdc;
    mu_i = pi/18;
    convergance = 0;
    
    while(convergance == 0)
        J = sqrt(3)/2/L*sin(mu_i);
        Idc_avg = 1/2/L*(-sqrt(3)*cos(mu_i)+sqrt(3));
        
        M = Idc_ref-Idc_avg;
        if abs(M) < 0.01
            break;
        end
        mu_i = mu_i+inv(J)*M;
    end
    
    if mu_i < 0 || mu_i > pi/6
        mu_i = pi/12;
    end
    mu_prev = mu_i*ones(12,1);
    gamma_prev = pi/6*ones(12,1);
    phaseAng_prev = angle(Valphabeta(2*h+3)+1i*Valphabeta(2*h+4));
else
    mu_prev = mu_in;
    gamma_prev = gamma_in;
    phaseAng_prev = phaseAng_in;
end
%************

Rdc = Rdc_in;
center = 2*h+1;
convergance = 0;

for i = 1:2:h
    V = Valphabeta(center-i*2:center+i*2+1);
    while(convergance == 0)      
        [Ia Ib Ic fundCurrents phaseAng gamma mu Pload error] =...
            solveDiode12Model(V,L,R,Rdc,Ldc,i,Vdc,P,mu_prev,gamma_prev,phaseAng_prev);
        
        ratio = Pload/P;
        if abs(1-ratio)*100 > 1.01
            c = clock;
            Rdc = Rdc*ratio;
            fprintf(errorFile,'%2.0fh:%2.0fm:%2.0fs Changing Rdc from %1.2f to %1.2f\n',c(4),c(5),c(6),Rdc/ratio,Rdc);
        elseif error == 1
            Rdc=Rdc*1.1;
        else
            convergance = 1;
        end
        
        if error == 0
            mu_prev = mu;
            gamma_prev = gamma;
            phaseAng_prev = phaseAng;
        end
        
    end
    mu_prev = mu;
    gamma_prev = gamma;
    phaseAng_prev = phaseAng;
    convergance = 0;
end

fclose(errorFile);