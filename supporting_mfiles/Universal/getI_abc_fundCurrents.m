function [Ia Ib Ic fundCurrents] = getI_abc_fundCurrents(I,h)
global SymMtx

Ia = zeros(h,3);
Ib = zeros(h,3);
Ic = zeros(h,3);
fundCurrents = zeros(3,2);

count = 1;

for i = 1:h
    
    outputCurrent = SymMtx*[0;I(2*h+2+count)+1i*I(2*h+2+count+1);...
        I(2*h-count)+1i*I(2*h-count+1)];
    
    if i == 1
        outa_ref = abs(outputCurrent(1));
        outb_ref = abs(outputCurrent(2));
        outc_ref = abs(outputCurrent(3));
        Ia(1,1:end) = [1,100.0,0.0];
        Ib(1,1:end) = [1,100.0,0.0];
        Ic(1,1:end) = [1,100.0,0.0];
        for j = 1:3
            fundCurrents(j,1) = abs(outputCurrent(j));
            fundCurrents(j,2) = angle(outputCurrent(j))*180/pi;
        end
    else %could combine the Ia,Ib,and Ic vectors into an array if time permits to clean up the code...not sure if it makes it any faster...somehow doubt it.      
        Ia(i,1:end) = [i,abs(outputCurrent(1))/outa_ref*100, angle(outputCurrent(1))*180/pi];
        Ib(i,1:end) = [i,abs(outputCurrent(2))/outb_ref*100, angle(outputCurrent(2))*180/pi];
        Ic(i,1:end) = [i,abs(outputCurrent(3))/outc_ref*100, angle(outputCurrent(3))*180/pi];
    end
    
    count = count + 2;
end
