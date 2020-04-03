%this function takes in the vector referencing phase a, b, and c and then
%combines the harmonic spectra that contains the magnitude and phase with
%respect to the fundamental
function Vec_sv = getSVMag(Vec_a,Vec_b,Vec_c,numHarmonics)
global invSymMtx

count = 0;
Vec_sv = zeros(2*numHarmonics+1,1);

for i = 1:numHarmonics
    Va_phasor = Vec_a(i,1)*exp(1i*Vec_a(i,2)*pi/180);
    Vb_phasor = Vec_b(i,1)*exp(1i*Vec_b(i,2)*pi/180);
    Vc_phasor = Vec_c(i,1)*exp(1i*Vec_c(i,2)*pi/180);
    posnegSeq = invSymMtx*[Va_phasor;Vb_phasor;Vc_phasor];
    Vec_sv(numHarmonics+1+i) = abs(posnegSeq(2));
    Vec_sv(numHarmonics+1-i) = abs(posnegSeq(3));
    count = count + 2;
end