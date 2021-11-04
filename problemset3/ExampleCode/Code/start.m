clc;clear;
dbstop if error

%% Set all options

% Things to do
options.solvecL     = 'N';      % Solve only c and L
options.solveeq     = 'Y';      % Solve equilibrium
options.polfun      = 'N';
options.sim         = 'N';      % Solve simulation
options.solveIRF    = 'Y';
options.solveKS     = 'N';      % Solve Krussel-Smith

% Model options
options.AR1         = 'Y';      % Approx continuous AR1 process, if 'N' then Rouwenhurst discretization
options.ACdown      = 'Y';      % If 'Y' then adjustment costs also for downwards movement
options.AC          = 'Y';      % If 'N' then adjustment costs are zero

% Tolerances, iterations
options.Nbell       = 2;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L

% Print / plot
options.print       = 'N';      % Print out c-solution convergence
options.plotSD      = 'N';      % Plot stationary distribution while solving equilibrium

%% Statespace parameters
glob.n          = [4,4];        % Number of nodes in each dimension
glob.nf         = [100,100];    % Number of points for k and z in histogram L
glob.curv       = 0.4;          % Curvature for k (1 is no curvature)
glob.spliorder  = [3,3];        % Order of splines (always use linear if productivity is discrete (not AR1))
glob.kmin       = 0.001;        % Lower bound on capital
glob.kmax       = 20.0;         % Upper bound on capital
glob.pzlb       = 0.005;        % Lower bound on probability of z
glob.Ne1        = 30;            % # of approx nodes of AR(1) iid shock in Expectation
glob.plb        = 0.001;        % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2        = 200;          % # of approx nodes of AR(1) iid shock in Approx of Q

% NOTE: Ne1 and Ne2 being very large only makes 'setup.m' run for longer,
% and not very much longer at that. They don't increase the dimensionality
% of *any* matrix

if glob.spliorder(2)>1 && strcmp(options.AR1,'N');fprintf('Error 1\n');return;end

%% Model parameters
% A. Outside parameters
glob.beta       = 0.99;
glob.rhoz       = 0.95;
glob.sige       = 0.05;
glob.A          = 1.00;

% B. Inside parameters (those you might want to change in a calibration
% exercise)
param.psi       = 2.4;
param.alpha     = 0.256;
param.nu        = 0.640;
param.delta     = 0.069;
param.eta       = 0.5;

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup(param,glob,options);
fprintf('Setup complete\n');

%% Solve only c and L for a given price p
switch options.solvecL
    case 'Y'
        p                   = 1.9;  % Conjectured value of p
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_cL(p,param,glob,options);
        fprintf('pin = %1.2f,\tpout = %1.2f\n',p,eq.p);
end

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolp        = 0.0001;             % Tolerance on price
        options.plb         = 0.5;              % Price lower bound
        options.pub         = 10;               % Price upper boud
        options.itermaxp    = 30;               % Max iterations of bisection
        options.eqplot      = 'N';
        options.eqprint     = 'Y';
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of p use old c as starting guess
        options.plotSD      = 'Y';              % If Y plot steady state distribution
        eq                  = solve_eq(param,glob,options);
end

save TEMP

%% Policy functions
load TEMP
switch options.polfun
    case 'Y'
        %__________________________________________________________________
        % 1. Varying Z for 3 levels of K
        kbar    = eq.L'*glob.sf(:,1);
        kplot   = [.5*kbar,kbar,1.25*kbar]';
        zplot   = glob.zgridf;
        kpolZ   = [];
        ipolZ   = [];
        for i = 1:numel(kplot);
            s           = gridmake(kplot(i),zplot);
            glob.Phi_Z  = splibas(glob.zgrid0,0,glob.spliorder(2),s(:,2));
            v           = solve_valfunc(eq.c,s,1,param,glob,options,1);
            kpolZ(:,i)  = v.Kp;
            ipolZ(:,i)  = v.I;
        end
        %__________________________________________________________________
        % 2. Varying K for 3 levels of Z
        zbar    = eq.L'*glob.sf(:,2);
        zplot   = [.5*zbar,zbar,1.25*zbar]';
        kplot   = glob.kgridf;
        kpolK   = [];
        ipolK   = [];
        for i = 1:numel(zplot);
            s           = gridmake(kplot,zplot(i));
            glob.Phi_Z  = splibas(glob.zgrid0,0,glob.spliorder(2),s(:,2));
            v           = solve_valfunc(eq.c,s,1,param,glob,options,1);
            kpolK(:,i)  = v.Kp;
            ipolK(:,i)  = v.I;
        end
        %__________________________________________________________________
        figure(round(1000*rand));
        subplot(2,2,1);plot(glob.kgridf,kpolK);title('1A. Capital');xlabel('Capital');legend('Low Z','Med Z','High Z');grid on;
        subplot(2,2,2);plot(glob.kgridf,ipolK);title('1B. Investment');xlabel('Capital');legend('Low Z','Med Z','High Z');grid on;
        subplot(2,2,3);plot(glob.zgridf,kpolZ);title('2A. Capital');xlabel('Productivity');legend('Low K','Med K','High K');grid on;
        subplot(2,2,4);plot(glob.zgridf,ipolZ);title('2B. Investment');xlabel('Productivity');legend('Low K','Med K','High K');grid on;
end

% Interesting: Constprod + s0

%% IRF
% Note: An IRF is only possible from the following (k,z)
% z must be the unconditional mean of the productivity process
% k must be the long-run capital associated with that z under no-shocks
load TEMP
% options.solveIRF = 'Y';
switch options.solveIRF
    case 'Y'
        if ~strcmp(options.AR1,'Y'),return,end;
        %__________________________________________________________________
        Ez                  = glob.Psszf'*glob.zgridf;
        %__________________________________________________________________
        % 1. Solve for long-run capital for z=E[z]=1 under no shocks
        options.simsolve    = 'N';
        options.constprod   = 'Y';
        options.IRF         = 'N';
        options.simplot     = 'N';
        T                   = 200;
        s0                  = gridmake(median(glob.kgridf),1);
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options);
        Kss                 = sim.Kt(end,:)';
        %__________________________________________________________________
        % 2. IRF - Start at [Kss,1] and solve with AR(1) shock in period 1
        options.simsolve    = 'Y';
        options.constprod   = 'N';
        options.IRF         = 'Y';
        options.simplot     = 'Y';
        T                   = 100;
        s0                  = [Kss,1];
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options);
end

%% Simulate model
switch options.sim
    case 'Y'
        %__________________________________________________________________
        % Initial vector of firms (can try different things)
        % 1. All k's average productivity
        s01                 = repmat(gridmake(glob.kgrid,glob.zgrid(round(glob.Nz/2))),500,1);
        % 2. Small k, all productivities - Growth
        s02                 = repmat(gridmake(glob.kgridf(10),glob.zgrid),1,1);
        % 3. Small k, all productivities - Growth
%         s04                 = repmat(gridmake(Kss,1*exp(-glob.sige)),5000,1);
        s04                 = repmat(gridmake(min(glob.kgrid),min(glob.zgridf)),10000,1);
        % 3. One individual point
        kbar                = eq.L'*glob.sf(:,1);
        zbar                = eq.L'*glob.sf(:,2);
        s03                 = [kbar,zbar];
        % Choose s0
        s0                  = s04;
        %__________________________________________________________________
        % Number of periods
        T                   = 400;
        % Options
        options.simsolve    = 'N';          % Solve problem each period (EXACT)
                                            % If N, then linearly
                                            % interpolate policies (FAST) -
        options.constprod   = 'N';          % Keep productivity constant
        options.simplot     = 'Y';
        options.IRF         = 'N';          % Fix as 'N'
        %__________________________________________________________________
        % Solve
        options.AC          = 'Y';
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options);
end




