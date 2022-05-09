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
% ECO384G Problem Set 6 (3, Spring 2022; 1, Chris), 1
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
%% Part 1: Write your own code to solve the problem using value function iteration
% Initialising intertemporal elasticity of substitution
sigma = 0.75;
% Initialising discount factor
beta = 0.95;
R = 1/beta - 0.00215;

% Setting iteration parameters
diff = 1;
epsi = 1e-4;

% Creating discrete search grid for state variables assets a and income Y
Amin = -20;
Amax = 60;
Na = 3*80;
Ymin = 2;
Ymax = 6;
Ny = 3;

% Discretisising state space
A = linspace(Amin, Amax, Na)';
Y = linspace(Ymin, Ymax, Ny)';

% Creating transition probability matrix
Ptran = [1/Ny 1/Ny 1/Ny];

% Initialising first guesses of value function
V = zeros(Na, Ny);

% Initialising first guesses of other variables
Vnew = zeros(Na, Ny);
Pindex = zeros(Na, Ny);

% Creating variable to keep track of iterations
iter = 0;

% Performing VFI
tic
while diff > epsi
  iter = iter + 1;
  for i = 1:Na
    for j = 1:Ny        
      C = A(i) + Y(j) - (A/R);
      % For values outside of the allowed grid, we set it to negative
      % infinity
      U = -inf(Na, 1);
      U(C > 0) = (1/(1 - 1/sigma))*(C(C > 0).^(1 - 1/sigma));
      Vtemp = U + beta*V*Ptran';
      [M, I] = max(Vtemp);
      Vnew(i, j) = M;
      Pindex(i, j) = I;
    end
  end  
  diff = max(max(abs(V - Vnew)));
  disp([iter diff]);
  
  % Saving current version of value function
  V = Vnew;
end
toc

% Plotting the value function
plot(A, V(:, 1), '-', 'Linewidth', 2);
hold on;
plot(A, V(:, 2), ':', 'Linewidth', 2);
plot(A, V(:, 3), '--', 'Linewidth', 2);
title('Value function');
xlabel('Current asset level');
ylabel('Value');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_1a_plot.png');
close(gcf);

% Plotting the policy function for assets a'
plot(A, A(Pindex(:, 1)), '-', 'Linewidth', 2);
hold on;
plot(A, A(Pindex(:, 2)), ':', 'Linewidth', 2);
plot(A, A(Pindex(:, 3)), '--', 'Linewidth', 2);
title('Policy function for next period assets');
xlabel('Current asset level');
ylabel('Next period`s asset level');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_1b_plot.png');
close(gcf);

% Plotting the policy function for consumption c'
Cs = repmat(A, 1, Ny) + repmat(Y', Na, 1) - [(A(Pindex(:, 1))/R) (A(Pindex(:, 2))/R) (A(Pindex(:, 3))/R)];
plot(A, Cs(:, 1), '-', 'Linewidth', 2);
hold on;
plot(A, Cs(:, 2), ':', 'Linewidth', 2);
plot(A, Cs(:, 3), '--', 'Linewidth', 2);
title('Policy for next period consumption');
xlabel('Current asset level');
ylabel('Next period`s consumption level');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'SouthEast');
hold off;
saveas(gcf, 'path\to\graphics\1_1c_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2: Speeding up VFI performed in part 1 using interpolation
% Initialising intertemporal elasticity of substitution
sigma = 0.75;
% Initialising discount factor
beta = 0.954;
R = 1/beta - 0.00215;

% Setting iteration parameters
diff = 1;
epsi = 1e-4;
itmax = 1000;
solverOptions = optimset('TolX', 1e-8);

% Creating utility function
global Utility;
Utility = @(c) (c > 0).*((1/(1 - 1/sigma))*c.^(1 - 1/sigma)) + (c <= 0).*-1e10;

% Creating discrete search grid for state variables assets a and income Y
Amin = -20;
Amax = 60;
Na = 1.5*30;
Ymin = 2;
Ymax = 6;
Ny = 3;

% Discretisising state space
A = linspace(Amin, Amax, Na)';
Y = linspace(Ymin, Ymax, Ny)';

% Creating transition probability matrix
Ptran = [1/Ny 1/Ny 1/Ny]';

% Initialising first guesses of value function
V = ones(Na, Ny);
ap_pol = zeros(Na, Ny);

% Initialising first guesses of other variables
U = ones(Na, 1);
Vnew = ones(Na, Ny);

% Creating variable to keep track of iterations
iter = 0;

% Performing VFI
tic
while diff > epsi && iter <= itmax
iter = iter + 1;
V = Vnew;
vFuture = beta*V*Ptran;
for i_a = 1:Na
  for i_y = 1:Ny      
    a = A(i_a);
    y = Y(i_y);
    ap_ub = a + y;
    [ap_pol(i_a,i_y), Vnew(i_a,i_y)] = fminbnd(@(ap) ...
      -value_of_ap(ap, i_a, i_y, Utility, A, Y, R, vFuture), A(1), ap_ub, solverOptions);
  end
end
Vnew = -Vnew;
diff = max(max(abs(V - Vnew)));
disp([iter diff]);
end
V = Vnew;
toc

% Creating axis showing asset levels found in VFI via interpolation
Adet = linspace(A(1), A(end), 200);

% Plotting value function
plot(Adet, interp1(A, V(:, 1), Adet, 'pchip', 'extrap'), '-', 'Linewidth', 2);
hold on;
plot(Adet, interp1(A, V(:, 2), Adet, 'pchip', 'extrap'), ':', 'Linewidth', 2);
plot(Adet, interp1(A, V(:,3), Adet, 'pchip', 'extrap'), '--', 'Linewidth', 2);
title('Value function')
xlabel('Current asset level')
ylabel('Value')
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location','Best');
hold off;
saveas(gcf, 'path\to\graphics\1_2a_plot.png');
close(gcf);

% Plotting the policy function for assets a'
plot(Adet, interp1(A, ap_pol(:, 1), Adet, 'pchip', 'extrap'), '-', 'Linewidth', 2);
hold on;
plot(Adet, interp1(A, ap_pol(:, 2), Adet, 'pchip', 'extrap'), ':', 'Linewidth', 2);
plot(Adet, interp1(A, ap_pol(:, 3), Adet, 'pchip', 'extrap'), '--', 'Linewidth', 2);
title('Policy function for next period assets');
xlabel('Current asset level');
ylabel('Next period`s asset level');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_2b_plot.png');
close(gcf);

% Computing consumption function
C = zeros(Na, Ny);
for i_a = 1:Na
  for i_y = 1:Ny      
    C(i_a, i_y) = Y(i_y) + A(i_a) - ap_pol(i_a, i_y)/R;
  end
end

% Plotting the policy function for consumption c'
plot(Adet, interp1(A, C(:, 1), Adet, 'pchip', 'extrap'), '-', 'Linewidth', 2);
hold on;
plot(Adet, interp1(A, C(:, 2), Adet, 'pchip', 'extrap'), ':', 'Linewidth', 2);
plot(Adet, interp1(A, C(:, 3), Adet, 'pchip', 'extrap'), '--', 'Linewidth', 2);
title('Policy funciton for next period consumption');
xlabel('Current asset level');
ylabel('Next period`s consumption level');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_2c_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 3: Assume income follows AR(1) process. Using Tauchen's procedure to discretise the AR(1) process and then performing VFI.
% Initialising intertemporal elasticity of substitution
sigma = 0.75;
% Initialising discount factor
beta = 0.954;
R = 1/beta - 0.00215;

% Setting iteration parameters
epsi = 1e-4;

% Creating discrete search grid for state variables assets a and income Y
Amin = -20;
Amax = 60;
Na = 3*80;
Ymin = 2;
Ymax = 6;
Ny = 3;

% Discretisising state space
A = linspace(Amin, Amax, Na)';
Y = linspace(Ymin, Ymax, Ny)';

% Setting various values for income AR(1) process that Tauchen's procedure
% will use
rho = [0.9 0.6];
sd = [0.5 1];
ss_val = 4;

% Initialising first guesses
Vs = zeros(Na, Ny, length(rho)*length(sd));
Cs = zeros(Na, Ny, length(rho)*length(sd));
As = zeros(Na, Ny, length(rho)*length(sd));
Ptran = zeros(Ny, Ny, length(rho)*length(sd));
P = zeros(Na*Ny, Ny, length(rho)*length(sd));

Agrid = repmat(A, 1, Ny);
AAgrid = reshape(Agrid, Na*Ny, 1);
Ygrid = repmat(Y, 1, Na)';
YYgrid = reshape(Ygrid, Na*Ny, 1);

CCgrid = AAgrid + YYgrid - A'./R;
UUgrid = -inf(Na*Ny, Na);
UUgrid(CCgrid > 0) = (CCgrid(CCgrid > 0).^(1-1/sigma))/(1 - (1/sigma));

% Keeping track of iterations
iter = 0;

% Performing Tauchen's discretisation process then VFI via vectorisation
% inside
tic
for i = 1:length(rho)
  for j = 1:length(sd)
    diff = 1;
    iter = iter + 1;
    disp([iter rho(i) sd(j)]);
    % Using Tauchen's procedure to discretise AR(1) process of income
    Ptran(:, :, iter) = tauchen1986(rho(i), sd(j), ss_val,Y');
    P(:, :, iter) = [repmat(Ptran(1, :, iter), Na, 1); repmat(Ptran(2, :, iter), Na, 1);...
      repmat(Ptran(3, :, iter), Na, 1)];
    VVgrid = zeros(Na*Ny,1);
    while diff > epsi
      VVtemp = UUgrid + beta*P(:, :, iter)*reshape(VVgrid, Na, Ny)';
      [VVnew, I] = max(VVtemp, [], 2);
      diff = max(abs(VVgrid - VVnew));
      VVgrid = VVnew;   
    end
    Vs(:, :, iter) = reshape(VVgrid, Na, Ny);
    Pindex = reshape(I, Na, Ny);
    As(:, :, iter) = [A(Pindex(:, 1)) A(Pindex(:, 2)) A(Pindex(:, 3))];
    Cs(:, :, iter) = Agrid + Ygrid - As(:, :, iter)/R;        
  end
end
disp(iter);
toc

% Plotting value function
plot(A, Vs(:, 1, 1), '-', 'Linewidth', 2);
hold on;
plot(A, Vs(:, 2, 1), ':', 'Linewidth', 2);
plot(A, Vs(:, 3, 1), '--', 'Linewidth', 2);
title('Value function assuming $\rho = 0.9, \sigma_{\epsilon} = 0.5$');
xlabel('Current asset level');
ylabel('Value');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location','Best');
hold off;
saveas(gcf, 'path\to\graphics\1_3a_plot.png');
close(gcf);

% Plotting the policy function for assets a'
plot(A, As(:, 1, 1), '-', 'Linewidth', 2);
hold on;
plot(A, As(:, 2, 1), ':', 'Linewidth', 2);
plot(A, As(:, 3, 1), '--', 'Linewidth', 2);
title('Policy function assuming $\rho = 0.9, \sigma_{\epsilon} = 0.5$');
xlabel('Current asset level');
ylabel('Next period`s asset level');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_3b_plot.png');
close(gcf);

% Plotting the policy function for consumption c'   
plot(A, Cs(:, 1, 1), '-', 'Linewidth', 2);
hold on;
plot(A, Cs(:, 2, 1), ':' , 'Linewidth', 2);
plot(A, Cs(:, 3, 1), '--', 'Linewidth', 2);
title('Policy function assuming $\rho = 0.9, \sigma_{\epsilon} = 0.5$');
xlabel('Current asset level');
ylabel('Next period`s consumption level');
legend('Y = 2', 'Y = 4', 'Y = 6', 'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_3c_plot.png');
close(gcf);

% Plotting policy function for consumption c' assuming different values of \rho and SD
plot(A, Cs(:, 1, 1), '-', 'Linewidth', 2);
hold on;
plot(A, Cs(:, 1, 2), ':', 'Linewidth', 2);
plot(A, Cs(:, 1, 3), '--', 'Linewidth', 2);
plot(A, Cs(:, 1, 4), '-.', 'Linewidth', 2);
title('Consumption policy function at $y = 4$, varying $\rho, \sigma_{\epsilon}$');
xlabel('Current asset level');
ylabel('Next period`s consumption level');
legend('$\rho = 0.9, \sigma_{\epsilon}=0.5$;',...
  '$\rho = 0.9, \sigma_{\epsilon} = 1$;',...
  '$\rho=0.6, \sigma_{\epsilon} = 0.5$;',...
  '$\rho = 0.6, \sigma_{\epsilon} = 1$;',...
  'Location', 'Best');
hold off;
saveas(gcf, 'path\to\graphics\1_3d_plot.png');
close(gcf);
%==========================================================================