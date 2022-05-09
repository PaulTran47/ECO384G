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

% This function calculates the utility for value function iteration via
% interpolation.
%==========================================================================

%==========================================================================
%=====
% NOTE
%=====
% Utility function should be declared in your primary programme.
%=========
% END NOTE
%=========
function val = value_of_ap(ap, i_a, i_y, Utility, A, Y, R, vFuture)
  val = Utility(Y(i_y) + A(i_a) - ap/R) + interp1(A, vFuture, ap, 'pchip');
end
%==========================================================================