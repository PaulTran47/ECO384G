%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% Comments that are important will be between the sub-section label:
%=====
% NOTE
%=====
% Important note here
%=========
% END NOTE
%=========

% This function implements Tauchen's 1986 procedure to discretise an AR(1)
% process. As in the paper, it is assumed that the grid is equidistant.

% Chris Boehm
% January 2014
%==========================================================================

%==========================================================================
%=====
% NOTE
%=====
% Arguments include:
% rho: Persistence of AR(1) process
% sig: Standard deviation of the INNOVATION to the process
% ss_val: mean of the process (steady state value)
% grid: grid the process should be discretised on
%=========
% END NOTE
%=========

function P = tauchen1986(rho, sig, ss_val, grid)  
  N = length(grid);
  step = grid(2) - grid(1);
  P = zeros(N, N);
  for j=1:N
    P(j, 1) = normcdf((grid(1) + step/2 - (1 - rho)*ss_val - rho*grid(j))/sig, 0, 1);
    for k = 2:(N - 1)
      P(j, k) = normcdf((grid(k) + step/2 - (1 - rho)*ss_val - rho*grid(j))/sig, 0, 1)...
        - normcdf((grid(k) - step/2 - (1 - rho)*ss_val - rho*grid(j))/sig, 0, 1);
    end
    P(j, N) = 1 - normcdf((grid(N) - step/2 - (1 - rho)*ss_val - rho*grid(j))/sig, 0, 1);
  end
end