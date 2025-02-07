function [idx, V, C] = SLC(X, Y, conf, P0)

[N, D] = size(Y);

if N==0
    idx=[]; V=[]; C = []; return
end

% Construct adjacency and laplace

M = 20;
[A, L] = Laplace_new(X, conf.delta, M);
M = size(L,1);

% Initialization
V = zeros(N,D); 
iter = 1;  
tecr = 1; 
C = zeros(M,D); 
E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D); 
gamma = conf.gamma;

%%
while (iter < conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) % 

    % E-step. 
    E_old = E;
    [P, E] = get_P(Y,V, sigma2, gamma, conf.a);
    if iter==1 && ~isempty(P0)
        P=P0;
    end  

    E = E+conf.lambda*trace(C'*L*C); 
    tecr = abs((E-E_old)/E);

    % M-step. 
    P = max(P, conf.minP);
    C = (A'.*repmat(P', [M, 1])*A + 2*conf.lambda*sigma2*L)\(A'.*repmat(P', [M, 1])*Y);

    % Update V and sigma^2
    V = A*C;
    Sp = sum(P);
    sigma2 = sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    % Update gamma
    numcorr = length(find(P > conf.theta));
    gamma = numcorr/size(X, 1);
    
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter = iter + 1;
end

%%
idx = find(P > conf.theta);

%%
function [P, E]=get_P(Y, V, sigma2 ,gamma, a)
% GET_P estimates the posterior probability and part of the energy.
D = size(Y, 2);
temp1 = exp(-sum((Y-V).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P = temp1./(temp1+temp2);
E = P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;
