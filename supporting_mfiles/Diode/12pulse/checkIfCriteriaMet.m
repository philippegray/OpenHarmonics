function [pass Vcurr] = checkIfCriteriaMet(MonitorList,Vprev,tol)
global numHarmonics OpenDSSFile

lenHarm = 2*numHarmonics+1;
exportMonitors(MonitorList);
Vcurr = zeros(lenHarm,length(MonitorList)/3);
pass = 1;

for i = 1:3:length(MonitorList)
    
    %does char(MonitorList(i)) because the objects for the abc phases are
    %all in sequence in the MonitorList structure.
    Va = csvread([OpenDSSFile,'_Mon_',char(MonitorList(i)),'.csv'],1,2,...
        [1,2,numHarmonics,3]);
    Vb = csvread([OpenDSSFile,'_Mon_',char(MonitorList(i+1)),'.csv'],1,2,...
        [1,2,numHarmonics,3]);
    Vc = csvread([OpenDSSFile,'_Mon_',char(MonitorList(i+2)),'.csv'],1,2,...
        [1,2,numHarmonics,3]);

    %Function combines the phase abc harmonic specta vectors into the space
    %vector format from -ve to +ve sequence harmonic h, that I have been
    %acustomed to using up to now
    Vcurr(1:end,(i-1)/3+1) = getSVMag(Va,Vb,Vc,numHarmonics);
end

compareMtx = Vcurr-Vprev;
for i = 1:2*numHarmonics+1
    for j = 1:length(MonitorList)/3
        if Vcurr(i,j) >= 1 %need this if statement for the even harmonics that are close to 0 but can vary from 10^-13 to 10^-12 for example each iteration which obviously doesn't matter which value...we don't care and it's a junk value.
            if abs(compareMtx(i,j)/Vprev(i,j))*100 > tol 
                pass = 0;
                break;
            end
        end
    end
    if pass == 0
        break;
    end
end

