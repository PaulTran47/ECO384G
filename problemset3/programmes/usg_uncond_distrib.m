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
% this function computes the unconditional distribution P associated with
% transition probability matrix TPM.
% END NOTE
%=========

%==========================================================================
function P = uncond_distrib(TPM,tol)
% P = uncond_distrib(TPM)
% P = uncond_distrib(TPM,tol) accepts a tolerance for the precision of P (default 1e-8). 

if nargin<2
  tol = 1e-8;
end

K = size(TPM,1);
dist_tpm = 1;
P = ones(1,K)/K;
while dist_tpm>tol
  Pnew = P*TPM;
  dist_tpm = max(abs(P(:)-Pnew(:)));
  P = Pnew;
end
P = P(:);