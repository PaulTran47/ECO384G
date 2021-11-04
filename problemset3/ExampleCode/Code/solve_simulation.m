function sim = solve_simulation(s0,T,eq,val,param,glob,options)

%% UNPACK
zgrid       = glob.zgrid;
c           = eq.c;
fspace      = glob.fspace;
Nz          = glob.Nz;
zub         = max(zgrid);
zlb         = min(zgrid);
rho         = glob.rhoz;
sigma       = glob.sige; 

%% SIMULATION

% Initial state
st          = s0;

% Storage (firms)
ni          = size(st,1); 
Kt          = zeros(T,ni);
Zt          = zeros(T,ni);  
Yt          = zeros(T-1,ni);
Nt          = zeros(T-1,ni);
Dt          = zeros(T-1,ni);
v1t         = zeros(T-1,ni);
It          = zeros(T-1,ni);
Ct          = zeros(T-1,ni);

% Storage (aggregate)
Kat         = zeros(T,1);

% Initial productivities
switch options.AR1
    case 'N'
        for i = 1:size(st,1);
            st(i,2) = find(glob.zgrid==st(i,2));    % Convert into indexes
        end
        Zt(1,:)     = glob.zgrid(st(:,2)');  
    case 'Y'
        switch options.IRF 
            case 'Y'
                Zt(1,:)     = exp(sigma)*st(:,2)'; 
                st(:,2)     = Zt(1,:)';
            case 'N'
                Zt(1,:)     = st(:,2)';
        end     
end 

% Initial capital
Kt(1,:)     = st(:,1)';
Kat(1,:)    = sum(Kt(1,:)'); 

% Policy interpolant
cKp         = funfitxy(fspace,glob.sf,eq.v.Kp); 
cY          = funfitxy(fspace,glob.sf,eq.v.Y);
cN          = funfitxy(fspace,glob.sf,eq.v.N);
cD          = funfitxy(fspace,glob.sf,eq.v.D);
cI          = funfitxy(fspace,glob.sf,eq.v.I);
cC          = funfitxy(fspace,glob.sf,eq.v.C);
cv1         = funfitxy(fspace,glob.sf,eq.v.v1);

% Draw productivity shocks
rng(999);
PS          = rand(ni,T); 

% Option: Constant productivity simulation
if strcmp(options.constprod,'Y');
    P       = eye(glob.Nz,glob.Nz);
else 
    P       = glob.P;        
end

% Simulation
for t = 2:T;
    % 1. Update state
    stm1                = st;
    % 2. Pull out state variables
     K               = stm1(:,1);
    switch options.AR1
        case 'Y'
            Z               = stm1(:,2);
        case 'N'
            Z               = glob.zgrid(stm1(:,2));
    end
    % 2. Compute values and policies given states stm1
    switch options.simsolve
        case 'Y'    
            % Recompute basis matrix for continuation values
            glob.Phi_Z      = splibas(glob.zgrid0,0,glob.spliorder(2),Z);
            % Solve problem (include 1 at end so don't compute E(v))
            v               = solve_valfunc(c,[K,Z],1,param,glob,options,1);
        case 'N'
            X               = funeval([cKp,cY,cN,cD,cI,cC,cv1],fspace,[K,Z]);
            v.Kp            = max(min(X(:,1),max(glob.kgrid)),min(glob.kgrid));
            v.Y             = X(:,2);
            v.N             = X(:,3);
            v.D             = X(:,4);
            v.I             = X(:,5);
            v.C             = X(:,6);
            v.v1            = X(:,7);
    end
    % 3. Evolve Z stochastically
    switch options.AR1
        case 'Y'
            switch options.constprod
                case 'Y'
                    % Keep constant
                    Zp          = Z;
                case 'N'
                    switch options.IRF
                        case 'N'
                            % Draw continuously distributed shock
                            Zp  = exp(rho*log(Z) + sigma*randn(size(Z)));
                        case 'Y'
                            % Evolve without shock
                            Zp  = exp(rho*log(Z));
                    end
            end
            Zp                  = max(min(Zp,zub),zlb);
            st                  = [v.Kp,Zp];
            Zt(t,:)             = st(:,2)';
        case 'N'
            Ptrans              = P(st(:,2),:);             % Transition probabilities for each individual
            Pcumsum             = cumsum(Ptrans')';
            x                   = PS(:,t);                  % Draw productivity shocks
            x                   = repmat(x,1,Nz);                  
            PPP                 = (Pcumsum>=x);  
            PPP                 = cumsum(PPP')';            % Numbers elements greater than x
            [~,J]               = find(PPP==1);
            Jp                  = J;
            st                  = [v.Kp,Jp];
            Zt(t,:)             = zgrid(st(:,2)');
    end
    %______________________________________________________________________
    % Record *next period state* for K
    Kt(t,:)             = st(:,1)';
    %______________________________________________________________________
    % Record period t objects, these depend on...
    Yt(t-1,:)           = v.Y';
    Nt(t-1,:)           = v.N';
    Dt(t-1,:)           = v.D';
    v1t(t-1,:)          = v.v1';
    It(t-1,:)           = v.I';
    Ct(t-1,:)           = v.C';
    %______________________________________________________________________
    % Aggregates
    Kat(t)              = sum(Kt(t,:)');
    Yat(t)              = sum(Yt(t-1,:)');
    Nat(t)              = sum(Nt(t-1,:)');
    Dat(t)              = sum(Dt(t-1,:)');
    Iat(t)              = sum(It(t-1,:)');
    Cat(t)              = sum(Ct(t-1,:)');
    
    fprintf('t: %3i\tAgg K: %1.3f\n',t,Kat(t)); 
end

%% PACK UP OUTPUT OBJECT {sim}
sim.Kt      = Kt;
sim.Zt      = Zt;
sim.Nt      = Nt;
sim.Dt      = Dt;
sim.v1t     = v1t;
sim.Kat     = Kat;

%% Figures
if strcmp(options.simplot,'Y')

    varlist     = {'Kt','Zt','Nt','It'};
    titlelist   = { 'A. Capital','B. Productivity','C. Labor','D. Investment'}; 

%     h4 = figure(round(rand*1000));
%     tvec = 0:1:T-1;
%     for x = 1:4;
%         subplot(2,2,x);
%         eval(sprintf('X=%s;',char(varlist(x))));
%         if x<=2
%             plot(tvec,X,'-','LineWidth',2);grid on;hold on; 
%         else
%             plot(tvec(1:end-1),X,'-','LineWidth',2);grid on;hold on; 
%         end
%         xlabel('Periods','fontsize',16);
%         title(titlelist(x),'fontsize',16);
%         set(gca,'fontsize',16);
%         xlim([0,T-1]);
%     end
    
    varlist     = {'Kat','Cat','Nat','Iat'};
    titlelist   = { 'A. Capital','B. Consumption','C. Labor','D. Investment'}; 
    
    h5 = figure(round(rand*1000));
    tvec = 0:1:T-1;
    for x = 1:4;
        subplot(2,2,x);
        eval(sprintf('X=%s;',char(varlist(x))));
        if x<=1
            plot(tvec,X,'-','LineWidth',2);grid on;hold on; 
        else
            plot(tvec(1:end-1),X(2:end),'-','LineWidth',2);grid on;hold on;  
        end
        xlabel('Periods','fontsize',16);
        title(titlelist(x),'fontsize',16);
        set(gca,'fontsize',16);
        xlim([0,T-1]);
    end
    
end

