function [v,jac] = solve_valfunc(c,s,P,param,glob,options,xxx)
c1 = c(1:end/2);
c2 = c(end/2+1:end);
%__________________________________________________________________________
% Solve problem 
B                       = menufun('bounds',s,[],P,param,glob,options); 
obj                     = @(kp)valfunc(c2,s,kp,P,param,glob,options);
Kp                      = goldenx(obj,B(:,1),B(:,2));
[v1,N,Y,I,D,AC,C,Phi_KpZ] = valfunc(c2,s,Kp,P,param,glob,options);
%__________________________________________________________________________
% Compute v2 and jacobian if requested
v2  = [];
if (nargin<=6)  % When simulating etc, don't need this, nargin is useful
    v2          = glob.Emat*c1;
end
if (nargout==2)
    jac         = [ glob.Phi,    -glob.beta*Phi_KpZ;
                    -glob.Emat,              glob.Phi ];
end
%__________________________________________________________________________
% Packup output
v.v1    = v1;
v.v2    = v2;
v.Kp    = Kp;
v.N     = N;
v.Y     = Y;
v.D     = D;
v.I     = I;
v.AC    = AC;
v.C     = C;




