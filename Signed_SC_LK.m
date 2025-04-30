% Signed SC-LK code

% Step 1- Load SSBM data
load 'SSBM_data.mat'

% Step 2- Run Signed SC-LK
k=5;
Alpha_range = unique([0.01, 0:0.05:1]);
rng(1); % For reproducibility
tic
for NL = 1:5
    Data_level = eval(sprintf('SSBM_%d', NL));    
    for is=1:size(Data_level,1)
        D=squeeze(Data_level(is,:,:,:));
        SSBM=permute(D, [2, 3, 1]);
        for ar=1:size(Alpha_range,2)
            Alpha=Alpha_range(ar);
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
                [U,V]=eig(Lbns_sys);
                [eigenvalue, index] = sort(diag(V));
                 k_smallest_pos_eigenvalues = eigenvalue(1:k); 
                 Ui = U(:, index(1:k));
                 Ui_norm = Ui./sqrt(sum(Ui.^2, 2));
                 Ui_norm_M(:,:,i)= Ui_norm;
                 UiUit_M(:,:,i)=Ui_norm*Ui_norm';
            end
            L_sum=sum(L_M,3);
            Ui_sum=sum(UiUit_M,3);
            L_mod=L_sum-(Alpha*(Ui_sum));
            [U_mod,V_mod]=eig(L_mod);
            [eigenvalue_mod, index_mod] = sort(diag(V_mod));
            k_smallest_eigenvalues_mod = eigenvalue_mod(1:k); 
            U_modk = U_mod(:, index_mod(1:k));
            U_modk_norm = U_modk./sqrt(sum(U_modk.^2, 2));
            [idx, C] = kmeans(U_modk_norm,k, 'Distance', 'sqEuclidean');
            NMI_values(ar)=getNMI(idx,GT);
            ARI_values(ar)=rand_index(idx,GT,'adjusted');
        end
        [~, y]=max(NMI_values);
        Optimal_Alpha(is,NL)=Alpha_range(y);
        NMI_value(is,NL)=NMI_values(y);
        ARI_value(is,NL)=ARI_values(y);
    end
    sprintf('SSBM_%d is done', NL)
end
toc