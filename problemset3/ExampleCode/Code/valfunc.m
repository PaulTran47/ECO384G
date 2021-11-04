function [v1,N,Y,I,F,AC,C,Phi_KpZ] = valfunc(c2,s,Kp,P,param,glob,options)
%__________________________________________________________________________
% Compute flow payoff
K           = s(:,1);
Z           = s(:,2);
F           = menufun('F',[K,Z],Kp,P,param,glob,options);
%__________________________________________________________________________
% Create basis matrices for continuation value
Phi_Kp      = splibas(glob.kgrid0,0,glob.spliorder(1),Kp);
Phi_KpZ     = dprod(glob.Phi_Z,Phi_Kp);    
%__________________________________________________________________________
% Compute value
v1          = F + glob.beta*Phi_KpZ*c2;       
%__________________________________________________________________________
% If requested, provide other policies
if nargout>1
    N       = menufun('labor',[K,Z],Kp,P,param,glob,options);
    Y       = menufun('output',[K,Z],Kp,P,param,glob,options);
    I       = menufun('investment',[K,Z],Kp,P,param,glob,options);
    AC      = menufun('costs',[K,Z],Kp,P,param,glob,options);
    C       = Y - I - AC;
end   