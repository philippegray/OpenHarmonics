function sv = getV_sv(filename_a,filename_b,filename_c,numHarmonics)

Va = csvread(filename_a,1,2,[1,2,numHarmonics,3]);
Vb = csvread(filename_b,1,2,[1,2,numHarmonics,3]);
Vc = csvread(filename_c,1,2,[1,2,numHarmonics,3]);

sv = getSVMag(Va,Vb,Vc,numHarmonics);