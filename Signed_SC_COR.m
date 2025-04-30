% Signed SC-COR code, use "Signed_Co_reg_func.m"

% Step 1- Load SSBM data
load 'SSBM_data.mat'

% Step 2- Run Signed SC-COR
epsilon=1e-4; k=5;
Alpha_range = unique([0.01, 0:0.05:1]);
rng(1); % For reproducibility
tic
for NL = 1:5
    Data_level = eval(sprintf('SSBM_%d', NL));    
    for is=1:size(Data_level,1)
        D=squeeze(Data_level(is,:,:,:));
        SSBM=permute(D, [2, 3, 1]);
        [N,~,M] = size(SSBM);
        for ar=1:size(Alpha_range,2)
            alpha=Alpha_range(ar);
            [NMI_V,ARI_V]=Signed_Co_reg_func(SSBM,alpha,k,epsilon,GT);
            NMI_values(ar)=NMI_V;
            ARI_values(ar)=ARI_V;
        end
        [~, y]=max(NMI_values);
        Optimal_Alpha(is,NL)=Alpha_range(y);
        NMI_value(is,NL)=NMI_values(y);
        ARI_value(is,NL)=ARI_values(y);
    end
    sprintf('SSBM_%d is done', NL)
end
toc