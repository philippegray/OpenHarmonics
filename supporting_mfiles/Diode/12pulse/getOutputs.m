function getOutputs(numHarmonics,MonitorList,sys,VSolve,ISolve) %I should have an "outputs" structure.
global DSSCircuit OpenDSSFile invSymMtx

%terms 1->numHarmonics = negative sequence magnitude
%terms numHarmonics+1 = dc magnitude
%terms numHarmonics+2 -> 2*numHarmonics+1 = postive sequence magnitude
%terms 2*numHarmonics+2 = dc-side dc-value.
%outputs = zeros(2*numHarmonics+1,length(MonitorList)/3*2);
outputs = struct();

%I think that 
%filename prefix
harmonicOrders = zeros(numHarmonics*2+1,1);

count = -numHarmonics;
for i = 1:numHarmonics*2+1
    harmonicOrders(i) = count;
    count = count + 1;
end

prefix = [OpenDSSFile,'_Mon_'];
Vmag = zeros(2*numHarmonics+1,1);
Imag = zeros(2*numHarmonics+1,1);
    
%currents injected from converter  - IEEE13Nodeckt_Mon_conv1a
for i = 1:sys.numConverters
    name_a = [prefix,sys.converterArr(i).name,'a.csv'];
    name_b = [prefix,sys.converterArr(i).name,'b.csv'];
    name_c = [prefix,sys.converterArr(i).name,'c.csv'];
    data_a = csvread(name_a,1,2,[1,2,numHarmonics,5]);
    data_b = csvread(name_b,1,2,[1,2,numHarmonics,5]);    
    data_c = csvread(name_c,1,2,[1,2,numHarmonics,5]);    
    
    newField = ['V_',sys.converterArr(i).name];
    outputs.(newField) = [harmonicOrders,getV_sv(name_a,name_b,name_c,numHarmonics)];
    
    Ia = data_a(1:end,3:4);
    Ib = data_b(1:end,3:4);
    Ic = data_c(1:end,3:4);
    
    newField = ['I_',sys.converterArr(i).name];
    outputs.(newField) = [harmonicOrders,getSVMag(Ia,Ib,Ic,numHarmonics)];

end

index = i;
%voltages of PCC corresponding to the converters
for i = 1:length(VSolve)
    filename_a = [prefix,'v_bus',char(VSolve(i)),'a.csv'];
    filename_b = [prefix,'v_bus',char(VSolve(i)),'b.csv'];
    filename_c = [prefix,'v_bus',char(VSolve(i)),'c.csv'];
    newField = ['V_bus',char(VSolve(i))];
    outputs.(newField) = [harmonicOrders,getV_sv(filename_a,filename_b,filename_c,numHarmonics)];
end

index = index+i;
%currents of lines
for i = 1:length(ISolve)
    
end

outputFile = 'outputs.mat';
save(outputFile,'outputs');