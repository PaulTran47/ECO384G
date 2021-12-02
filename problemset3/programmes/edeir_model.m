%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% Answers to question parts that don't involve code can be found at the
% bottom of the programme, in the section ``Questions asked in problemset x
% that don't involve code".

% Text answers to question parts that involve code will be between the
% sub-section label:
%=======
% ANSWER
%=======
% Answer here
%===========
% END ANSWER
%===========

% Comments that are important will be between the sub-section label:
%=====
% NOTE
%=====
% Important note here
%=========
% END NOTE
%=========
% ECO384G Problem Set 3, 3 and 4
% Paul Le Tran, plt377
% 1 December, 2021
%==========================================================================

%=====
% NOTE
%=====
% The following code is of Martin Uribe's design, which serves as the
% baseline for the numerical solutions to this problem set. Specifically,
% this function computes analytically the derivatives of the equilibrium
% conditions of the EDEIR open economy model presented in Chapter 4 of
% ``Open Economy Macroeconomics" by Martin Uribe and Stephanie
% Schmitt-Grohe.
% END NOTE
%=========

%==========================================================================
function edeir_model
% The equilibrium conditions take the form
% E_t f(yp,y,xp,x) =0, 
% where x is a state vector and y a control vector, and  'p' denotes
% next-period. The derivatives of this  system evaluated at the
% deterministic steady state  are denoted fx, fxp, fy, and fyp up to first
% order accuracy, the  dynamics of this system are governed by the system
% xp = hx * x +ETASHOCK epsp
% y = gx * x, 
% where ETASHOCK is a matrix of paramters and epsp~(0,I). 
% The variance-covariance matrix of the forcing term is denoted varshock=ETASHOCK*ETASHOCK';

% Output: edeir_model.m produces the file edeir_model_num_eval.m 
% containing symbolic expressions for f, fx, fxp, fy, fyp, and  ETASHOCK 

% Calls: anal_deriv.m and first_order_print2f.m

filename = 'edeir_model'; %output file name (it will have a suffix _num_eval in the .m file used for numerical evaluation)

sourcename = mfilename('fullpath'); %carries the name of the current file and its full path.

approx = 1; %Order of approximation desired

% Define parameters
syms RSTAR BETA DBAR DELTA ALPHA  PHI  RHO  STD_EPS_A ETATILDE SIGG OMEGA PSSI

% Define endogenous  variables 
variables = {'la', 'c', 'k',  'kfu',  'h', 'd', 'output',  'ivv', 'tb', 'tby' ,'ca','cay', 'r'};
% Create a symbol for the variable in period t and period t+1
aux1 =length(variables); 
for aux2=1:aux1; 
  eval(['syms ' variables{aux2} ' ' variables{aux2} 'p']);
end
% Exogenous Shocks
shocks = {'a'};
% Create symbols for the shocks in period t and period t+1
aux1 =length(shocks); 
for aux2=1:aux1; 
  eval(['syms ' shocks{aux2} ' ' shocks{aux2} 'p']);
end
% Standard deviations of exogenous shocks
Standard_Deviations  = {'ETATILDE'};

% Equilibrium conditions. The symbols e1, e2, ... denote equation 1, equation2, ...
% Evolution of debt
e1 = -dp + (1+r)*d + c + ivv + PHI/2*(kp -k)^2 - output; 

% Output
e2 = -output + a*k^ALPHA*h^(1-ALPHA);

% FOC w.r.t. consumption
e3 = -la + (c-h^OMEGA/OMEGA)^(-SIGG); 

% FOC w.r.t. h
e4 = -h^(OMEGA-1)+ (1-ALPHA) * a * (k/h)^ALPHA;

% FOC w.r.t. debt
e5 = -la + BETA * (1+rp) * lap;

% FOC w.r.t. capital
e6 = -la* (1+PHI*(kp-k)) + BETA * lap * (1-DELTA + ALPHA * ap * (kp/hp)^(ALPHA-1) + PHI * (kfup-kp));

% Country premium
e7 = -rp + RSTAR + PSSI * (exp(dp-DBAR) -1);

% Investment
e8 = -ivv +kp - (1-DELTA)*k;

% Trade balance
e9 = -tb + output - c - ivv;

% Trade-balance-to-output ratio
e10 = -tby + tb/output;

% Current account
e11 = -ca -dp + d;

% Current-account-to-output ratio
e12 = -cay + ca/output;

% Evolution of TFP
e13 = -log(ap) + RHO * log(a);

% Make kfu=kp
e14 = -kfu+kp;

% Create function f
f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14];

% Create the vector of state and control variables
statevar = [];
controlvar = [];

% Define the vector of controls, y, and states, x
statevar = [statevar d  r k  a];
statevarp = p_it(statevar);%the program p_it (a function at the bottom of this code, creates the state in period t+1

controlvar = [controlvar c ivv output h la kfu tb tby ca cay];
controlvarp = p_it(controlvar);

% Pick variables around which f will be  log-linearized (f will be linearized around all other variables)
log_linearize = [ k a c ivv output h la kfu ];
log_linearizep = p_it(log_linearize);

% Variables to substitute from levels to logs
log_linearize = [log_linearize log_linearizep];
f = subs(f, log_linearize, exp(log_linearize));

% Compute analytical derivatives of f
[fx,fxp,fy,fyp]=anal_deriv(f,statevar,controlvar,statevarp,controlvarp,approx);%the program anal_deriv.m is available at: http://www.columbia.edu/~mu2166/1st_order/1st_order.htm

% Make f and its derivatives a function of the level of its arguments rather than the log
f = subs(f, log_linearize, log(log_linearize));
fx = subs(fx, log_linearize, log(log_linearize));
fy = subs(fy, log_linearize, log(log_linearize));
fxp = subs(fxp, log_linearize, log(log_linearize));
fyp = subs(fyp, log_linearize, log(log_linearize));

% Eliminate future variables
cu = [statevar controlvar];
cup = [statevarp controlvarp];

f = subs(f, cup,cu,0);
fx = subs(fx, cup,cu,0);
fy = subs(fy, cup,cu,0);
fxp = subs(fxp, cup,cu,0);
fyp = subs(fyp, cup,cu,0);

STANDARD_DEVIATIONS = []; aux1 =length(shocks); 
for  aux2=1:aux1; 
STANDARD_DEVIATIONS = [STANDARD_DEVIATIONS sym(Standard_Deviations{aux2})]
end
nshocks = numel(shocks);
nstates = numel(statevar);
ETASHOCK(nstates,nshocks) = sym(0);
for aux2=1:nshocks
aux1 = find(statevar==shocks(aux2))
ETASHOCK(aux1,aux2) = STANDARD_DEVIATIONS(aux2);
end

varshock = ETASHOCK*ETASHOCK';
% Print derivatives to file <filename>_num_eval.m'  for model evaluation
first_order_print2f(filename,fx,fxp,fy,fyp,f,ETASHOCK,statevar,controlvar,sourcename); %this script can be downloaded from http://www.columbia.edu/~mu2166/1st_order/1st_order.htm

% Positions of variables in state and control vectors
nstate = length(statevar);
for i=1:nstate
eval(['n' char(statevar(i)) ' = ' num2str(i) ';']);
end;

ncontrol = length(controlvar);
for i=1:ncontrol
eval(['n' char(controlvar(i)) ' = ' num2str(i) ';']);
end;


eval(['save ' filename '.mat'])

function  statevarp = p_it(statevar);
aux1=length(statevar);
statevarp = [];
for aux2 = 1:aux1
x = [char(statevar(aux2)) 'p'];
statevarp = [statevarp  sym(x)];
end
%==========================================================================