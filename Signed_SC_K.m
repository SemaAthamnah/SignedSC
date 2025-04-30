% Signed SC-K code

% Step 1- Load SSBM data
load 'SSBM_data.mat'

% Step 2- Run Signed SC-K
k=5;
rng(1); % For reproducibility
tic
for NL = 1:5
    Data_level = eval(sprintf('SSBM_%d', NL));    
    for is=1:size(Data_level,1)
        D=squeeze(Data_level(is,:,:,:));
        SSBM=permute(D, [2, 3, 1]);
        M=size(SSBM,3);
        clear L_M Ui_norm_M UiUit_M
        for i=1:M
            A=SSBM(:,:,i);
            Ap=A.*(A>0);
            An=-A.*(A<0);
            Dp=diag(sum(Ap,2));
            Dn=diag(sum(An,2));
            Dpar=diag(sum(abs(A),2));
            Lbs= Dp-A;
            Lbns_sys= Dpar^(-0.5)*Lbs* Dpar^(-0.5);
            L_M(:,:,i)=Lbns_sys;
            [Ui,V]=eigs(Lbns_sys,k, 'SA');
            Ui_norm = Ui./sqrt(sum(Ui.^2, 2));
            Ui_norm_M(:,:,i)= Ui_norm;
            UiUit_M(:,:,i)=Ui_norm*Ui_norm';
        end
        SumUiUit=sum(UiUit_M,3);
        [Ustar, Estar] = eigs(SumUiUit, k, 'LA'); 
        Ustar=Ustar./sqrt(sum(Ustar.^2, 2));
        [idx, C] = kmeans(Ustar,k, 'Distance', 'sqEuclidean');
        NMI_value(is,NL)=getNMI(idx,GT);
        ARI_value(is,NL)=rand_index(idx,GT,'adjusted');
    end
    sprintf('SSBM_%d is done', NL)
end
toc