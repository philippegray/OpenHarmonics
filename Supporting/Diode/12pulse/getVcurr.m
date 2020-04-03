function Vcurr = getVcurr(len)
global numHarmonics OpenDSSFile

exportMonitors(MonitorList);

Vcurr = zeros(2*numHarmonics+1,len/3);

for i = 1:len/3
    
end
Va = csvread([OpenDSSFile,'_Mon_',vsc_obj.name,'a.csv'],1,2,...
                [1,2,numHarmonics,3]);
Vb = csvread([OpenDSSFile,'_Mon_',vsc_obj.name,'b.csv'],1,2,...
                [1,2,numHarmonics,3]);
Vc = csvread([OpenDSSFile,'_Mon_',vsc_obj.name,'c.csv'],1,2,...
                [1,2,numHarmonics,3]);

count = 0;
for i = 1:h
    Vav = Va(i,1)*exp(1i*Va(i,2)*pi/180);
    Vbv = Vb(i,1)*exp(1i*Vb(i,2)*pi/180);
    Vcv = Vc(i,1)*exp(1i*Vc(i,2)*pi/180);
    pnArr = invSymMtx*[Vav;Vbv;Vcv];
    V(2*h+3+count:2*h+4+count) = [real(pnArr(2));imag(pnArr(2))];
    V(2*h-1-count:2*h-count) = [real(pnArr(3));imag(pnArr(3))];
    count = count + 2;
end
V = V/Vbase;

