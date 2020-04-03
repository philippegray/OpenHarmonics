function harmVec = harmonicsVector(numHarmonics)

harmVec = '';
for i = 1:numHarmonics
   harmVec = [harmVec,[' ',num2str(i)]];
end

harmVec = harmVec(2:end);

