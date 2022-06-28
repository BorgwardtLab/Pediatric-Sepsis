addpath('/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/')
addpath('/cluster/work/borgw/JK/virtual_environments/spss_cluster/Utilities/');

Input = load('/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/data/phy_dem_jk.mat');
X = double(Input.Dat);
z = double(Input.Cov);
for k = 2:6
    disp(k)
    
    [~, ~, U] = KCC(X, z, k, k);
    
    csvwrite(strcat('/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/data/phy_dem_jk_',string(k),'.csv'), U)
end
