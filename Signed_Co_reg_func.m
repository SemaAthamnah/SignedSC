function [NMI_V,ARI_V]=Signed_Co_reg_func(SSBM,alpha,k,epsilon,GT)

% Signed SC-COR function 

[N,~,M] = size(SSBM);
Ui = zeros(N,k,M);
clear L_M Ui
for i=1:M
    A=SSBM(:,:,i);
    Ap=A.*(A>0);
    An=-A.*(A<0);
    Dp=diag(sum(Ap,2));
    Dn=diag(sum(An,2));
    Dpar=diag(sum(abs(A),2));

    Lbs= Dp-A;
    L_M(:,:,i)=Dpar^(-0.5)*Lbs* Dpar^(-0.5);
    [U Vi] = eigs(L_M(:,:,i),k,'SA'); % Min eigenvalues (smallestreal)
    Ui(:,:,i) = U./sqrt(sum(U.^2, 2)); 

    %the objective function as the sum of the diagonal elements 
    % of the eigenvalue matrix (sum(diag(Vi))) is based on 
    % Rayleigh quotient properties in spectral clustering and Laplacian-based optimization.
    %objval(i,1) = sum(diag(Vi));
end

prev_obj = inf; 
iter = 1;
EV=[];
N1=[];
N2=[];
clear prev_Ustar prev_Ui curr_Ustar curr_Ui
prev_Ustar = inf(300,5);
prev_Ui=Ui;
while true
    %fprintf('Iteration %d\n', iter);
    % Fix Ui, solve for U*
    SumUiUit = zeros(N, N);
    for j = 1:M
        SumUiUit = SumUiUit + (prev_Ui(:, :, j) * prev_Ui(:, :, j)');
    end
    [Ustar, Estar] = eigs(SumUiUit, k, 'LA'); % Max eigenvalues (largestreal)
    curr_Ustar=Ustar./sqrt(sum(Ustar.^2, 2));
    UstarUstart = curr_Ustar * curr_Ustar';

    % Fix U*, solve for Ui
    for j = 1:M
        %[Ui(:, :, j), E] = eigs(L_M(:, :, j) - (alpha * UsUst), k, 'SA'); % Min eigenvalues
        [U, E] = eigs(L_M(:, :, j) - (alpha * UstarUstart), k, 'SA'); 
        curr_Ui(:, :, j)=U./sqrt(sum(U.^2, 2));
        objval(j, iter + 1) = sum(diag(E));
    end

    % Check convergence (N1,N2)-relative Frobenius norm difference 
    %N1(iter)=norm(curr_Ustar-prev_Ustar,'fro')/norm(prev_Ustar,'fro');
    %N2(iter)=norm(curr_Ui-prev_Ui,'fro')/norm(prev_Ui,'fro');
    % 
    % if (N1(iter)< epsilon) && (N2(iter)< epsilon)
    %     disp('Convergence reached')
    %     break;
    % end
    prev_Ustar = curr_Ustar;
    prev_Ui = curr_Ui;
    % iter = iter + 1;

    % Check convergence (Objval)
    curr_obj = sum(objval(:, iter + 1));
    EV(iter)=curr_obj;

    if abs(prev_obj - curr_obj) < epsilon
        break;
    end
    prev_obj = curr_obj;
    iter = iter + 1;
end

rng(1); % For reproducibility
[idx, ~] = kmeans(curr_Ustar, k, 'Distance', 'sqEuclidean');
NMI_V=getNMI(idx,GT);
ARI_V=rand_index(idx,GT,'adjusted');

end