function Zi_arr = rotMtx(lenHarm,h,gamma,phase,Zi)

RotationMatrix = zeros(3*lenHarm+1,3*lenHarm+1,13);
Zi_arr = zeros(3*lenHarm+1,13);

gamma_sum = 0;

for i = 1:13
    count = 1;
    for j  = -h:1:h
        RotationMatrix(count:count+1,count:count+1,i) = ...
            [cos(j*(gamma_sum-phase)),-sin(j*(gamma_sum-phase));...
            sin(j*(gamma_sum-phase)),cos(j*(gamma_sum-phase))];
        count = count + 2;
    end
    for j = 0:h
        RotationMatrix(count:count+1,count:count+1,i) = ...
            [cos(j*(gamma_sum-phase)),-sin(j*(gamma_sum-phase));...
            sin(j*(gamma_sum-phase)),cos(j*(gamma_sum-phase))];
        count = count + 2;
    end
    if i < 13
        gamma_sum = gamma_sum + gamma(i);
    end
    Zi_arr(1:end,i) = RotationMatrix(1:end,1:end,i)*Zi;
end
