addpath([pwd '/Utilities/']);
data_path = 'data';
Input = load(strcat(data_path, '/Cov_AgeSexEth_NormalImputation.mat'));
X_physio = double(Input.Physio);
X_contextual = double(Input.Contextual);
X_clinical = double(Input.Clinical);
X_proteome = double(Input.Proteome);
z = double(Input.Cov);

for k = 19:20
    
    [~, ~, U_physio] = KCC(X_physio, z, k, k);
    csvwrite(strcat(data_path, '/KCC_Cov_AgeSexEth_physio_NormalImputation_K',string(k),'.csv'), U_physio)
    [~, ~, U_contextual] = KCC(X_contextual, z, k, k);
    csvwrite(strcat(data_path, '/KCC_Cov_AgeSexEth_contextual_NormalImputation_K',string(k),'.csv'), U_contextual)
    [~, ~, U_clinical] = KCC(X_clinical, z, k, k);
    csvwrite(strcat(data_path, '/KCC_Cov_AgeSexEth_clinical_NormalImputation_K',string(k),'.csv'), U_clinical)
    [~, ~, U_proteome] = KCC(X_proteome, z, k, k);
    csvwrite(strcat(data_path, '/KCC_Cov_AgeSexEth_proteome_NormalImputation_K',string(k),'.csv'), U_proteome)
end

