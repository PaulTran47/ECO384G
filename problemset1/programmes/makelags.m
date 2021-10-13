%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% This function returns the original series plus lags.
%==========================================================================

%==========================================================================
function x=makelags(x0,k)
[n1,n2]=size(x0);
if n2>1
    x0=x0';
end

x(:,1)=x0(1+k:length(x0));
for j=1:k
    x=[x x0(1+k-j:length(x0)-j)];
end
return
%==========================================================================