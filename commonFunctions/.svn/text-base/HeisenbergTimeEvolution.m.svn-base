function [ spin, spin_ex ] = HeisenbergTimeEvolution( H , dt , x0 , B , steps , err )
% HeisenbergTimeEvolution: 
%% Time Evolution
L = length(x0);
N = numel(B(1,:));
%exact solution

[V,E] = eig(full(H));

%exakt time evoulution operator
T_ex =  V * diag(exp(-1i * dt * diag(E))) * V';

%% time steps
spin = zeros(N,steps); % Spin konfigurationen lanczos
spin_ex = zeros(N,steps); % -||- exakt

x0_ex = x0; %exakter zustandsvektor

%keyboard
for k = 1:steps
    
    %exakt
    x0_ex = T_ex * x0_ex; % zeitentwicklung exakt
    x0_ex1 = repmat(x0_ex,1,N);
    
    spin_ex(:,k) = sum(B .* (x0_ex1 .* conj(x0_ex1))); % spin update exakt
    
    
    %lanczos time evolution
    [xn, conv_flag]  = lanczosTimeEvolution_myversion( H , L , err , L , x0 , dt);
    
    
    if conv_flag == 1
        
        xn1 = repmat(xn,1,N);
        spin(:,k) = sum(B .* (xn1.*conj(xn1)));
        
        x0=xn;
    end
  
end


