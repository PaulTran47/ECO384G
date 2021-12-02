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
% this function computes a number of first- and second-moments from
% transition probability matrix TPM amd matrix of states x. TPM is of order
% n-by-n and x of order n-by-nx. Each column of x is understood as a
% different random variable, all with transition probability matrix TPM.
% END NOTE
%=========

%==========================================================================
function [Ex,Vx,SDx,Corrx,Scorrx,Covx,Acovx,P] = moments_tpm(TPM,x,T,P);
% [Ex,Vx,SDx,Corrx,Scorrx] = moments_tpm(TPM,x,P)
% Ex is the unconditional mean of x (order 1-by-nx)
% Vx is the unconditional variance of x (order 1-by-nx)
% SDx is the unconditional standard deviation of x (order 1-by-nx)
% Corrx is the unconditional correlation of x (order nx-by-nx)
% and Scorrx  are the autocorrelations of x of orders 1 to T (order T-by-nx), default T=1;
% Inputs: 
% TPM is the transition probability matrix (order n-by-n)
% x is a matrix of nx random variables (order n-by-nx)
% [Ex,Vx,SDx,Corrx,Scorrx] = moments_tpm(TPM,x,T) allows to set  the maximum order of autocorrelations (default 1)
% [Ex,Vx,SDx,Corrx,Scorrx] = moments_tpm(TPM,x,T,P) allows to give the  unconditional distribution  implied by TPM, which is necessary to compute unconditional first and second moments.  If not provided, the programs computes it. 

% Note: The same output is produced by the program mt.m using a nonvectorized (loop) method to compute second  moments. 

% Number of states
n = size(TPM,1);

% Number of random variables
nx = size(x,2);

if nargin<4
  %Unconditional distribution 
  dist_tpm = 1;
  P = ones(n,1)/n;
  TTPM = TPM';
  while dist_tpm>1e-8
    Pnew = TTPM*P;
    dist_tpm = max(abs(P(:)-Pnew(:)));
    P = Pnew;
  end
end
clear TTPM Pnew

if nargin<3
  T=1;
end

% Means
Ex = P'*x;

% Deviations from means
% dx = x-repmat(Ex,n,1);
dx = bsxfun(@minus, x, Ex); 

% Variances
Vx = P'*(dx.^2);

% Standard deviations
SDx = Vx.^(1/2);

% Variacne-Covariance matrix
bsxfun(@times,dx,P);
Covx = ans'*dx;

% Correlations
Corrx = Covx./(SDx'*SDx);

% Autocovariances of order t=1 to T, Ex_t*x_t+j'
% Serial correlation of order t=1 to T
aux = TPM;
for j=1:T
  bsxfun(@times,dx,P);
  Acovx(1:nx,1:nx,j) = ans'*aux*dx;
  Scorrx(j,1:nx) = diag(Acovx(:,:,j))'./Vx;
  aux = aux*TPM;
end