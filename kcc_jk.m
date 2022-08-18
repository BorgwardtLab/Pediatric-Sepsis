addpath('/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/')
addpath('/cluster/work/borgw/JK/virtual_environments/spss_cluster/Utilities/');

f0 = '/cluster/work/borgw/SPSS/MultiOmicsAnalysis/ConsensusClustering/data/'
f2 = '_dem_jk'
fs = ["phy","con","cli","pro"]

for i = 1:length(fs)

    f1 = fs(i)
    fprintf(f1)

    Input = load(strcat(f0,f1,f2,'.mat'));
    X = double(Input.Dat);
    z = double(Input.Cov);

    for k = 2:6
        for d = k:6
    
            disp(k)
            [~, ~, U] = KCC(X, z, k, d);
            csvwrite(strcat(f0,f1,f2,'_k_',string(k),'_d_',string(d),'.csv'), U)

        end
    end

end