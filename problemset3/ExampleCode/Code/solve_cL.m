function eq  = solve_cL(P,param,glob,options)

%% A. Globals 
s           = glob.s;  
sf          = glob.sf;
kgrid       = glob.kgrid;
ns          = size(s,1);

%% B. Compute equilibrium objects that depend on p
% ----- None -----

%% Initialise guesses (if val.cresult has an old guess in it, use that)
c1old       = zeros(ns,1);
c2old       = zeros(ns,1);
c           = options.cresult;
if ~isempty(c) && strcmp(options.Loadc,'Y')
    c1old    = c(1:ns);
    c2old    = c(ns+1:end);
end
cold        = [c1old;c2old];
totaltic    = tic;

%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfunc(cold,s,P,param,glob,options); 
    % 2. Update c
    c1          = glob.Phi\full(v.v1); 
    c2          = glob.Phi\full(v.v2); 
    c           = [c1;c2];
    % 3. Compute distance and update
    dc          = norm(c-cold)/norm(cold); 
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
end

%% Newton iterations
if strcmp(options.print,'Y');
    fprintf('~~~~~ Newton iterations ~~~~~\n');
end
eq.flag.cconv = false;
for citer = (1:options.Nnewt)
    % 1. Compute values
    [v,jac]     = solve_valfunc(cold,s,P,param,glob,options);
    % 2. Update c 
    c1old       = cold(1:glob.Ns); 
    c2old       = cold(glob.Ns+1:end);
    c           = cold - jac\([glob.Phi*c1old - full(v.v1) ;
                               glob.Phi*c2old - full(v.v2)]);  
    % 3. Compute distances and update
    dc          = norm(c-cold)/norm(cold);
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
    % 4. Check convergence
    if (dc<options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv,break,end;
end

%% Solve again on a finder grid for k
glob.Phi_Z      = glob.Phi_Zf; 
v               = solve_valfunc(c,sf,P,param,glob,options,1);

%% Compute stationary distribution
Kp              = min(v.Kp,max(kgrid));
fspaceergk      = fundef({'spli',glob.kgridf,0,1});
QK              = funbas(fspaceergk,Kp); 
QZ              = glob.QZ;
Q               = dprod(QZ,QK); 

% [vv,dd]         = eigs(Q');
% dd              = diag(dd);
% Lv              = vv(:,dd==max(dd));
% L               = Lv/sum(Lv);
L               = ones(size(Q,1),1);
L               = L/sum(L);

for itL = (1:options.itermaxL);
    Lnew    = Q'*L;  
    dL      = norm(Lnew-L)/norm(L);  
    if (dL<options.tolL),break,end;
    if mod(itL,100)==0 
        if strcmp(options.print,'Y')
            fprintf('dL:\t%1.3e\n',dL);
        end
    end
    L       = Lnew;
end

% Plot stationary distribution
if strcmp(options.plotSD,'Y');
    H = figure(888);
    set(H,'Pos',[1          41        1920         964]);
    Jk  = numel(glob.kgridf);
    Jz  = numel(glob.zgridf);
    Lz  = kron(eye(Jz),ones(1,Jk))*L; 
    Lk  = kron(ones(1,Jz),eye(Jk))*L;
    % Marginal K
    subplot(2,2,1);
    plot(glob.kgridf,Lk,'o-');title('Stationary Capital Dist - Lk');grid on;
    % Marginal Z
    subplot(2,2,2);
    plot(exp(glob.zgridf),Lz,'o-');title('Stationary Prod Dist - Lz');grid on;
    eq.Lk   = Lk; 
    eq.Lz   = Lz;
    % Joint (K,Z) - Surface plot
    subplot(2,2,4)
    Lmat    = reshape(L,Jk,Jz);
    Zmat    = repmat(glob.zgridf',Jk,1); 
    Kmat    = repmat(glob.kgridf,1,Jz);
    if Lk(end)<0.001;
        Kub     = glob.kgridf(find(cumsum(Lk)>0.98,1,'first'));
    else
        Kub     = max(glob.kgridf);
    end
    Zub     = glob.zgridf(find(cumsum(Lz)>0.98,1,'first')); 
    mesh(Zmat,Kmat,Lmat,'LineWidth',2);
    xlabel('Productivity - z');ylabel('Capital - k');title('Joint Distribution');
    xlim([min(glob.zgridf),Zub]);
    ylim([min(glob.kgridf),Kub]); 
    zlim([0,max(max(Lmat))]);
    drawnow;
end

%% Compute aggregates and implied price p
Ya      = L'*v.Y;
Ia      = L'*v.I;
Ka      = L'*sf(:,1);
ACa     = L'*v.AC;
Ca      = Ya - Ia - ACa;     
p       = 1/Ca;

%% Pack-up output
eq.v    = v;
eq.c    = c;
eq.p    = p;
eq.L    = L;
eq.Ya   = Ya;
eq.Ca   = Ca;
eq.Ia   = Ia;
eq.ACa  = ACa;
eq.Ka   = Ka;
eq.Q    = Q;

end

