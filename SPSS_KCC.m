addPath('/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/')
addpath([pwd '/Utilities/']);

Input = load('data/ProteomeView_AgeSex.mat');
X = double(Input.Proteome);
z = double(Input.Cov);
for k = 2:6
    disp(k)
    
    [~, ~, U] = KCC(X, z, k, k);
    
    csvwrite(strcat('./data/TwoViewsKCC_proteome_K_',string(k),'.csv'), U)
end

