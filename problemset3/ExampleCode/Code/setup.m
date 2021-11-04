function [param,glob] = setup(param,glob,options)

%% State space for idiosyncratic productivity
% One persistent shock
Nz              = glob.n(2);
switch options.AR1
    case 'Y' 
        logzub              = norminv(1-glob.pzlb,0,glob.sige/sqrt(1-glob.rhoz^2));
        logzlb              = norminv(glob.pzlb,0,glob.sige/sqrt(1-glob.rhoz^2));
        zlb                 = exp(logzlb);
        zub                 = exp(logzub);
        zgrid               = nodeunif(Nz,zlb,zub);
        % For comparison compute Rouwenhurst grid and Pssz
        [P,zgridRouw,Pssz]  = setup_MarkovZ(Nz,glob.sige,glob.rhoz,1);
    case 'N'
        [P,zgrid,Pssz]      = setup_MarkovZ(Nz,glob.sige,glob.rhoz,1);         
end
zgrid0          = zgrid;

%% State space for endogenous variable k
Nk              = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
kgrid           = nodeunif(Nk,glob.kmin.^curv(1),glob.kmax.^curv(1)).^(1/curv(1));  % Adds curvature
kgrid0          = kgrid;    % Save for computing basis matrices in valfunc.m:line9

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',kgrid,0,glob.spliorder(1)},...
                         {'spli',zgrid,0,glob.spliorder(2)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
kgrid           = s(s(:,2)==s(1,2),1); 
zgrid           = s(s(:,1)==s(1,1),2);
Nk              = size(kgrid,1); 
Nz              = size(zgrid,1);

%% Compute expectations matrix
switch options.AR1
    case 'Y'
        Ne              = glob.Ne1;
        pvec            = nodeunif(Ne,glob.plb,1-glob.plb);     % Make an equi-spaced grid in probabilities
        e               = norminv(pvec,0,glob.sige);            
        w               = normpdf(e,0,glob.sige);               % Invert normal for shocks
        w               = w/sum(w);                             % Compute pdf of shocks
        iNe             = ones(Ne,1);                           
        iNs             = ones(Ns,1);
        gfun            = @(z,e) max(min(exp(glob.rhoz*log(z)+e),max(zgrid)),min(zgrid));   % Constrained to lie within nodes
        g               = gfun(kron(s(:,2),iNe),kron(iNs,e));
        Phi             = funbas(fspace,[kron(s(:,1),iNe),g]);
        Ikronw          = kron(eye(Ns),w');
        glob.Emat       = Ikronw*Phi;        
    case 'N'
        Phi             = funbas(fspace,s);
        glob.Emat       = kron(P,speye(Nk))*Phi;
end

%% Construct fine grid for histogram
kgridf          = nodeunif(glob.nf(1),glob.kmin.^glob.curv(1),glob.kmax.^glob.curv(1)).^(1/glob.curv(1));
Nkf             = size(kgridf,1);

switch options.AR1
    case 'Y'
        zgridf  = nodeunif(glob.nf(2),min(zgrid),max(zgrid));
    case 'N'
        zgridf  = zgrid;
end        
Nzf             = size(zgridf,1);
sf              = gridmake(kgridf,zgridf);
Nsf             = size(sf,1);

glob.kgridf     = kgridf;
glob.zgridf     = zgridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QZ matrix for approximation of stationary distribution
switch options.AR1
    case 'Y' 
        % 1. AR(1)
        Ne              = glob.Ne2;
        pvec            = nodeunif(Ne,glob.plb,1-glob.plb);         % Make an equi-spaced grid in probabilities
        e               = norminv(pvec,0,glob.sige);                % Invert normal for shocks
        w               = normpdf(e,0,glob.sige);                   % Compute pdf of shocks
        w               = w/sum(w);                                 % Normalise
        fspaceZ         = fundef({'spli',zgridf,0,1});              % Linear interpolant
        QZ              = zeros(Nsf,Nzf);
        P               = zeros(Nzf,Nzf);                           % P constructed so can compute steady state Psszf and compare to Pssz
        for i = 1:Ne;
            g           = gfun(sf(:,2),e(i));
            QZi         = funbas(fspaceZ,g);
            QZ          = QZ + w(i)*QZi;
            P           = P  + w(i)*funbas(fspaceZ,gfun(zgridf,e(i)));
        end
        glob.QZ         = QZ;
        % For plotting
        Psszf           = P^1000;
        Psszf           = Psszf(1,:)';
        glob.Psszf      = Psszf;
        if strcmp(options.plotSD,'Y')
            figure(round(1000*rand));
            subplot(1,2,1);
            plot(zgridRouw,Pssz,'bo-');grid on; hold on;title('A. Coarse Transition matrix Pssz - From Rouwenhurst'); 
            subplot(1,2,2);
            plot(zgridf,Psszf,'ro-');grid on;title('B. Fine stationary dist');
        end
    case 'N'
        % 2. Discrete productivity process
        glob.QZ         = kron(P,ones(Nkf,1)); 
end

%% Create one time only basis matrices
glob.Phi_Z      = splibas(zgrid0,0,spliorder(2),s(:,2));                % Used in Bellman / Newton computing expected values
glob.Phi_Zf     = splibas(zgrid0,0,spliorder(2),sf(:,2));               % Used when solving on fine grid
Phi_K           = splibas(kgrid0,0,spliorder(1),s(:,1));
glob.Phi        = dprod(glob.Phi_Z,Phi_K);                              % Used in Bellman / Newton updating of c

%% Declare additional global variables
glob.kgrid0     = kgrid0;
glob.kgrid      = kgrid;
glob.zgrid0     = zgrid0;
glob.zgrid      = zgrid;
glob.P          = P;
glob.Pssz       = Pssz;
glob.Nz         = Nz;
glob.Nk         = Nk;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;
%__________________________________________________________________________
end


        
        
        
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
   

