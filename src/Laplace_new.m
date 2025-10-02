function [W, L] = Laplace_new(X, delta, M)

corres = X;
[N,~] = size(corres);

rng(3);

[corres1, coridx] = unique(corres,'rows');
M = min(M, size(corres1,1));
con_ps = coridx(randperm(size(corres1,1),M));

W = repmat(corres,[1 1 M])-permute(repmat(corres(con_ps,:),[1 1 N]),[3 2 1]);
W = exp(-squeeze(sum(W.^2,2))./(delta^2)); % 

D = W(con_ps,:);
L = diag(sum(D,2)) - D;
WLW = D'*L*D; 

a=1;