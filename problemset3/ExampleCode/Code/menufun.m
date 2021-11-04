function out = menufun(flag,s,x,P,param,glob,options)
% Globals
    % None
% Parameters
nu          = param.nu;
psi         = param.psi;
alpha       = param.alpha; 
delta       = param.delta;
eta         = param.eta;
% Equilibrium variables and their dependents
w           = psi./P;

switch flag
    case 'bounds'
        kpmin   = ones(size(s,1),1)*min(glob.kgrid); 
        kpmax   = ones(size(s,1),1)*max(glob.kgrid);
        out     = [kpmin,kpmax];
    case 'output'
        K       = s(:,1);
        Z       = s(:,2);
        N       = Nfun(Z,K);
        out     = Yfun(Z,N,K); 
    case 'labor'
        K       = s(:,1);
        Z       = s(:,2);
        out     = Nfun(Z,K);
    case 'investment'
        K       = s(:,1);
        Kp      = x;
        out     = Ifun(K,Kp);
    case 'costs'
        K       = s(:,1);
        Kp      = x;
        out     = ACfun(K,Kp); 
    case 'F'
        K       = s(:,1);
        Z       = s(:,2);
        Kp      = x;
        N       = Nfun(Z,K);        % Labor demand
        Y       = Yfun(Z,N,K);      % Output
        Pi      = Y - w.*N;         % Operating profits
        I       = Ifun(K,Kp);       % Investment
        AC      = ACfun(K,Kp);      % Adjustment costs
        out     = P.*(Pi-I-AC);     % Final flow payoff
end

%__________________________________________________________________________
% NESTED FUNCTIONS
% 1. Labor demand from FOC(n')
function N = Nfun(Z,K)
    N = ((nu*Z.*K.^alpha)./w).^(1/(1-nu));
end        
% 2. Production function
function Y = Yfun(Z,N,K)
    Y = Z.*(K.^alpha).*(N.^nu);
end  
% 3. Production function
function AC = ACfun(K,Kp)
    switch options.AC
        case 'Y'
            switch options.ACdown
                case 'Y'
                    AC = eta*(Kp-K).^2; 
                case 'N'
                    AC = (Kp>K).*(eta*(Kp-K).^2); 
            end
        case 'N'
            AC = zeros(size(K));
    end
end
% 4. Investment
function I = Ifun(K,Kp)
    I = Kp - (1-delta)*K;
end
%__________________________________________________________________________


        
end
        
        
        
