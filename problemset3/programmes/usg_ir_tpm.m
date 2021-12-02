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
% this function computes impulse responses of the variables in x to a shock
% specified by s.
% END NOTE
%=========

%==========================================================================
function IR = usg_ir_tpm(PAI,uPAI,x,s,T)
% IR = ir_tpm(PAI,uPAI,x,s,T)
% PAI is the transition probability matrix associted with x
% uPAI is the unconditional probability  distribution of x
% x is a matrix with as many rows as the order of PAI and as many columns as the number of variables
% s is a vector of indices defining the information set in the first period. For example, if all that is known is that in the initial period the economy  is in states 2, 7, or 13, then s=[2 7 13];
 
if nargin<5
  T = 21; %number of periods
end

% Number of variables
v = size(x,2);

PAI0 = zeros(size(x(:,1)));
PAI0(s) = uPAI(s);
PAI0 = PAI0./sum(PAI0);

for t=1:T
  bsxfun(@(x,y) x.*y,x,PAI0);
  IR(t,1:v) = sum(ans);
  x = PAI*x;
end