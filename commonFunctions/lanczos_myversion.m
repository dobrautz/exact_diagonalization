function [ E,V3 ] = lanczos_myversion( H , varargin)
%lanczos_myversion Lanczos Algorithmus
% Input:  H     %Hamilton Operator
%         k     %welcher,wieviele Eigenwerte?
%         err   Genauigkeit
%         nmax  max. Anzahl lanczos schritte

% output:  E     Eigenwerte
%          V3    Eigenvektoren
%
%------------------------SVN Info------------------------------------------
% $Rev:: 61                             $: Revision der letzten Übertragung
% $Author:: dobrautz                    $: Autor der letzten Übertragung
% $Date:: 2013-06-11 18:24:48 +0200 (Di#$: Datum der letzten Übertragung
% -------------------------------------------------------------------------

% default values and struct call handling
if numel(varargin) == 0
    
    k = 1;
    err = 10^-5;
    nmax = 100;
    
elseif numel(varargin) == 1
    
    k = varargin{1}.nEigenvalues;
    err = varargin{1}.precision;
    nmax = varargin{1}.nSteps;
    
elseif numel(varargin) == 3;
    
    k = varargin{1};
    err = varargin{2};
    nmax = varargin{3};
    
end


if isempty(H) == 1
    
    E = [];
    
    V3 = [];
    
elseif sum(size(H))==2;
    
    E = H(1);
    V3 = 1;
    
elseif length(H) < 20
    
    
    if k > length(H)
        k = length(H);
    end
    
    if isreal(H)
        [V3,E] = eigs(H,k,'SA');
        E = diag(E);
    else
        [V3,E] = eigs(H,k,'SR');
        E=diag(E);
    end
    
else
    N = length(H);
    V3 = zeros(N,k);
    E = zeros(k,nmax);
    ortho_ind = zeros(1,nmax);
    
    if k > N
        error('k > N!')
    else
        
        %1. iteration: zufallsvektor
        xn = rand(N,1); %zufaellige koeffizienten
        
        xn = xn/norm(xn);
        En = real(xn'*H*xn);
        
        xn1 = H*xn - En*xn;
        
        
        V(:,1) = xn; T(1,1) = En;
        %E(1) = eigs(T,1,'SA');
        
        for j = 2:nmax
            
            % Recursion
            xn_1 = xn;
            xn = xn1;
            
            kn = norm(xn);
            
            xn = xn/norm(xn);
            En = real(xn'*H*xn);
            
            xn1 = H*xn - En*xn - kn*xn_1;
            
            T(j,j-1) = kn; T(j-1,j) = kn; T(j,j) = En;
            V(:,j)   = xn;
            
            % Reorthogonalize
            h1 = V' * xn1;
            
            if norm(h1) > 10*sqrt(eps)
                xn1 = xn1 - V * h1;
                xn1 = xn1 - V * (V' * xn1);
                ortho_ind(j) = 1;
            end
            
            
            if j > k
                E1 = eig(T);
                E(:,j) = min(E1);
                
                %
                if abs(E(k,j)-E(k,j-1))< err || kn < err
                    
                    break
                    
                end
            end
        end
        ortho_ind = sum(ortho_ind);
        
        %% Eigenwerte
        E = E(:,j);
        
        %% Eigenvektoren
        [V2,~] = eig(T);
        
        for f = 1:k
            V3(:,f) = V*V2(:,f); %grundzustandsvektor
            
        end
        
        
    end
    
end
