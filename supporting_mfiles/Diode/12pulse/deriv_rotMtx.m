function derivRotation = deriv_rotMtx(lenHarm,h,gamma,phase)

maxPos = 3*lenHarm+1;
derivRotation = zeros(maxPos,maxPos,13);
gamma_sum = 0;

for i = 1:13
    count = 1;
    for j  = -h:1:h
        derivRotation(count:count+1,count:count+1,i) = ...
            [-j*sin(j*(gamma_sum-phase)),-j*cos(j*(gamma_sum-phase));...
            j*cos(j*(gamma_sum-phase)) -j*sin(j*(gamma_sum-phase))];
        count = count + 2;
    end
    for j  = 0:1:h
        derivRotation(count:count+1,count:count+1,i) = ...
            [-j*sin(j*(gamma_sum-phase)),-j*cos(j*(gamma_sum-phase));...
            j*cos(j*(gamma_sum-phase)) -j*sin(j*(gamma_sum-phase))];
        count = count + 2;
    end
    if i < 13
        gamma_sum = gamma_sum + gamma(i);
    end
end
