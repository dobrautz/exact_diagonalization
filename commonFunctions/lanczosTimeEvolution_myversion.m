function [ xt_neu ,conv_flag] = lanczosTimeEvolution_myversion( H , k , t_err , nmax , xn , dt)
%lanczosTimeEvolution_myversion Summary of this function goes here
%   Detailed explanation goes here

N = length(H);

ortho_ind = zeros(1,nmax);

if k > N
    error('k > N!')    
else
    

xn = xn/norm(xn);
En = xn'*H*xn;
xn1 = H*xn - En*xn;


V(:,1) = xn; T(1,1) = En;

xt_alt = 0;
conv_flag = 0;

for j = 2:nmax
    
    % Recursion
    xn_1 = xn;
    xn = xn1;
       
    kn = norm(xn);
    
    xn = xn/norm(xn);
    En = xn'*H*xn;
    
    xn1 = H*xn - En*xn - kn*xn_1;
    
    T(j,j-1) = kn; T(j-1,j) = kn; T(j,j) = En;
    V(:,j)   = xn;
    
%     Reorthogonalize
%     h1 = V(:,1:j)'*xn1;
%     
%     if norm(h1) > 10*sqrt(eps)
%         xn1 = xn1 - V(:,1:j)*h1;
%         xn1 = xn1 - V(:,1:j)*(V(:,1:j)'*xn1);
%         ortho_ind(j) = 1;
%     end
    
    
    [G,E] = eig(T);
    
    %time evolution
    T_evo = G * diag(exp(-1i * dt * diag(E))) * G';
    
    xt_neu = ((V * T_evo) * V') * V(:,1);
    
    if all(abs(xt_alt - xt_neu) < t_err)
       conv_flag =1; 
       break
    end
    
    xt_alt = xt_neu;
    
end

end