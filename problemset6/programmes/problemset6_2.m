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
% ECO384G Problem Set 6 (3, Spring 2022; 1, Chris), 2 (Problem 1)
% Paul Le Tran, plt377
% 7 May, 2022
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';

% Setting text interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

cd(home_dir);
%==========================================================================

%==========================================================================
%% Iterating to obtain next period's durable goods level
%
% Initialising intertemporal elasticity of substitution
sigma = 0.75;
% Initialising discount factor
beta = 0.95;
r = 1/beta - 0.00215;
% Initialising depreciation rate
delta = 0.1;
% Initialising fixed cost
F = 3;
% Initialising price of durable good D
p = 1.5;

% Setting iteration parameters
epsi = 1e-4;

% Creating discrete search grid for state variables assets a
Amin = 0;
Amax = 1;
Na = 3*10;

% Setting income level
y = 5;

% Creating discrete search grid for durable goods D
Dmin = 0.01;
Dmax = 40;
Nd = 3*50;

% Creating discrete search grid for random depreciation stock e
emin = -0.05;
emax = 0.05;
Ne = 3;

% Discretisising state space
A = linspace(Amin, Amax, Na)';
D = linspace(Dmin, Dmax, Nd)';
e = linspace(emin, emax, Ne)';

% More discretisation
A1 = repmat(A', Ne, Nd);
a1 = A1(:);
e1 = repmat(e, Nd*Na, 1);
D1 = repmat(D', Ne*Na, 1);
d1 = D1(:);
ddgrid = D';
d2 = repmat(D', Na*Nd*Ne, 1);

% Creating transition probability matrix
Ptran = [1/Ne 1/Ne 1/Ne];

% Initialising array to store next period's D
a2 = zeros(Nd*Na*Ne, Nd);

% Performing iteration to obtain next period's D
for i = 1:size(a2, 1)
  for j = 1:size(a2, 2)
    if ddgrid(j) == d1(i)
      a2(i, j) = (y + a1(i) - p*(ddgrid(j) - d1(i)*(1 - delta - e1(i))))*(1 + r);
      d2(i, j) =ddgrid(j);
      if (y + a1(i) - p*(ddgrid(j) - d1(i)*(1 - delta-e1(i))))*(1 + r) <= Amin 
        d2(i, j) = (y + a1(i) - F - Amin/(1 + r))/p + d1(i)*(1 - delta - e1(i));
        a2(i, j) = Amin;
      else
        if (y + a1(i) - p*(ddgrid(j) - d1(i)*(1 - delta - e1(i))))*(1 + r) >= Amax
          d2(i, j) = (y + a1(i) - F - Amax/(1 + r))/p + d1(i)*(1 - delta - e1(i));
          a2(i, j) = Amax;
        end
      end
    else 
      a2(i, j) = (y + a1(i) - F - p*(ddgrid(j) - d1(i)*(1 - delta - e1(i))))*(1 + r);
      d2(i, j) = ddgrid(j);
      if (y + a1(i) - F - p*(ddgrid(j) - d1(i)*(1 - delta - e1(i))))*(1 + r) <= Amin 
        d2(i, j) = (y + a1(i) - F - Amin/(1 + r))/p + d1(i)*(1 - delta - e1(i));
        a2(i, j) = Amin;
      else
        if (y + a1(i) - F - p*(ddgrid(j) - d1(i)*(1 - delta - e1(i))))*(1 + r) >= Amax
          d2(i, j)=(y + a1(i) - F-Amax/(1 + r))/p + d1(i)*(1 - delta - e1(i));
          a2(i, j) = Amax;
        end
      end
    end           
  end
end

% Assuming Utility is of the form u(D) = ln(D)
U = log(d2);

max(d2,[], "all")
min(d2,[], "all")
%==========================================================================

%==========================================================================
%% Performing VFI via vectorisation
% Setting iteration parameters
diff = 1;

% Initialising first guesses of value and policy function
v = zeros(Ne, Na, Nd);
P = repmat(Ptran, Ne, 1);

tic
while diff > epsi
  Ev = P*v(:, :);
  Ev = Ev(:); 
  [vnew,b] = max(U + beta*Ev, [], 2);
  diff = max(abs(v(:) - vnew));
  v = reshape(vnew, Ne, Na, Nd);
end
toc

dp1 = zeros(length(b), 1);
ap1 = zeros(length(b), 1);

for i = 1:size(d2, 1)
  for j = 1:size(d2, 2)
    dp1(i) = d2(i, b(i));
    ap1(i) = a2(i, b(i));
  end
end

dp = reshape(dp1, Ne, Na, Nd);
ap = reshape(ap1, Ne, Na, Nd);
%==========================================================================

%==========================================================================
%% Plotting value and policy functions
plotvalue = find(A == max(A((A < median(A) + 0.2 & median(A) - 0.2))));

% Plotting value function
plot(D, reshape(v(1, plotvalue, :), 1, length(D)), '-', 'Linewidth', 2);
hold on;
plot(D, reshape(v(2, plotvalue, :), 1, length(D)), ':', 'Linewidth', 2);
plot(D, reshape(v(3, plotvalue, :), 1, length(D)), '--', 'Linewidth', 2);
title('Value function at asset level = 13.14');
xlabel('Current durable stock level');
ylabel('Value function');
legend('$\epsilon = -0.05$', '$\epsilon = 0$', '$\epsilon = 0.05$', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\2_1a_plot.png');
close(gcf);

% Plotting policy function for durable good D'
plot(D, reshape(dp(1, plotvalue, :), 1, length(D)), '-', 'Linewidth', 2);
hold on;
plot(D, reshape(dp(2, plotvalue, :), 1, length(D)), ':', 'Linewidth', 2);
plot(D, reshape(dp(3, plotvalue, :), 1, length(D)), '--', 'Linewidth', 2);
title('Policy function for durable good D at asset level 13.14');
xlabel('Current durable stock level');
ylabel('Next period`s durable stock level');
legend('$\epsilon = -0.05$', '$\epsilon = 0$', '$\epsilon = 0.05$', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\2_1b_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Performing simulation to find long-run/ergodic distribution of durable good D
% Setting total number of simulations to perform
T = 10000;

% Initialising random number generator
rng('default');
rng(1010);

% Creating matrix of random numbers
X = rand(T, 1);

% Creating variable for standard error
Se = zeros(T, 1);
Se(X >= 0 & X < 1/3) = -0.05;
Se(X >= 1/3 & X < 2/3) = 0;
Se(X >= 2/3 & X <= 1) = 0.05;

Se_index = zeros(T, 1);
Se_index(X >= 0 & X < 1/3) = 1;
Se_index(X >= 1/3 & X < 2/3) = 2;
Se_index(X >= 2/3 & X <= 1) = 3;

% Initialising durable goods and asset levels
Ds = zeros(T, 1);
As = zeros(T, 1);
bindex = reshape(b, Ne, Na, Nd);

% Assuming initial level (i.e., at period 0) of D is 0.01, A = e = 0
D0 = 0.01;
A0 = 0;
e0 = 0;
Ds(1) = D(bindex(2, 1, 1));

if Ds(1) == D0
  As(1) = (y + A0 - p*(Ds(1)-D0*(1 - delta - Se(1))))*(1 + r);
else
  As(1) = (y + A0 - F - p*(Ds(1) - D0*(1 - delta - Se(1))))*(1 + r);
end

% Performing simulation
for i = 2:T
  [~, ida] = min(abs(As(i - 1) - A));
  idd = find(D == Ds(i - 1));
  Ds(i) = D(bindex(Se_index(i), ida, idd));
  if Ds(i) == Ds(i - 1)
    As(i) = (y + As(i - 1) - p*(Ds(i) - Ds(i - 1)*(1 - delta - Se(i))))*(1 + r);
  else 
    As(i) = (y + As(i - 1) - F - p*(Ds(i) - Ds(i - 1)*(1 - delta - Se(i))))*(1 + r);
  end
end

% Plotting ergodic distribution of D through normalising histogram
Dh = histogram(Ds);
Dh.Normalization = 'probability';
title('Ergodic distribution of durable good D');
xlabel('Durable stock level');
ylabel('Frequency');
hold off;
saveas(gcf, 'path\to\graphics\2_1c_plot.png');
close(gcf);
%==========================================================================