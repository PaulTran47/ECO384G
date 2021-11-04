function [eq,val] = solve_eq(param,glob,options)

plb0    = options.plb;
pub0    = options.pub;
plb     = plb0;
pub     = pub0;

% Storage
pinvec          = zeros(options.itermaxp,1); 
poutvec         = zeros(options.itermaxp,1);
options.cresult = [];
tictic = tic;
for tt = (1:options.itermaxp)
    % 1. Update p 
    p               = (1/2)*(plb+pub);
    % 2. Solve economy given p
    eq              = solve_cL(p,param,glob,options);  
    options.cresult = eq.c;         % Save to use as starting guess
    % 3. Record output and print
    pinvec(tt)      = p;
    poutvec(tt)     = eq.p;
    if strcmp(options.eqprint,'Y') 
        fprintf('%2i. pin:\t%2.6f\tpout:\t%2.6f\tt:%2.1f\n',tt,p,eq.p,toc(tictic));
    end
    % 4. Set all flags
    d               = pinvec-poutvec;
    eq.flag.exist   = ~all(sign(d)==max(sign(d))); 
    eq.flag.equi    = (abs(pinvec(tt)-poutvec(tt))<options.tolp);
    eq.flag.down    = (pinvec(tt)>poutvec(tt));
    eq.flag.up      = (pinvec(tt)<poutvec(tt));
    % 5. Shift bounds
    plb             = (eq.flag.up)*p    + (eq.flag.down)*plb;
    pub             = (eq.flag.up)*pub  + (eq.flag.down)*p;
    % 6. Break if equilibrium
    if eq.flag.equi,break,end
    % 7. Plot option
    if strcmp(options.eqplot,'Y') 
       figure(888);       
       subplot(2,2,3);
       plot(pinvec(pinvec~=0),'bo','markersize',6,'markerfacecolor','b');hold on;grid on;
       plot(poutvec(poutvec~=0),'ro','markersize',6,'markerfacecolor','r');
       xlim([0,max(tt,12)]);
       legend('pin','pout','Location','NorthEast');
       title('Equilibrium');
       drawnow;  
    end
end
