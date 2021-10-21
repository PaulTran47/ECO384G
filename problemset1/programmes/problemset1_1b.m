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

% ECO384G Problem Set 1, 1b
% Paul Le Tran, plt377
% 1 September, 2021
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';
data_dir = 'path\to\data';

cd(home_dir);
%==========================================================================

%==========================================================================
%% Loading in data
cd(data_dir);
data1 = xlsread(append(data_dir, '\INDPRO.xls'));
data2 = xlsread(append(data_dir, '\RomerandRomerDataAppendix.xls'), 'DATA BY MONTH');
%==========================================================================

%==========================================================================
%% Creating time series
% Setting up the data to be the same data range as Romer^2 (2004) data.
% Following Romer^2 in that time series start at 1966m1.

% THIS IS IF YOU ARE USING FRED'S IP DATA.
ip_sa = data1(564:936, 1); % Included 1965m12 to have 1966m1 in diffs.
chg_ln_ip_sa = diff(log(ip_sa));
mp_shock = data2(1:end, 1);

% Making regressand follow sample period of 1970m1 - 1996m12.
chg_ln_ip_sa = chg_ln_ip_sa(49:end, :);

% Deleting datasets as no longer needed.
clear data1 data2;
%==========================================================================

%==========================================================================
%% Making lags for mp_shock
cd(home_dir);

% Matrix include current and lagged series
% This uses's Coibion's makelags.m
h = 48;
mp_shockl = makelags(mp_shock, h);

% Removing current vectors from lag matrix.
mp_shockl = mp_shockl(:, 2:end);
%==========================================================================

%==========================================================================
%% Creating simple regression
regressors = [ones(size(chg_ln_ip_sa)) mp_shockl];
b = regress(chg_ln_ip_sa, regressors);

% b_1 = coefficients of mp_shockl
b_1 = b(2:end, :);
%==========================================================================

%==========================================================================
%% Creating IRF
% For a moving average representation, the estimators are the IRFs.
irf = b_1;

cirf = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf(i) = cirf(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf = vertcat(a, cirf);
clear a;
%==========================================================================

%==========================================================================
%% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_ip_sa - y_hat;

% Creates a matrix of estimators from bootstrapping process
b_b_matrix = zeros(h + 1, B);
for i = 1:B
  resid_b = datasample(resid, length(resid)); % Reshuffling residuals
  y_b = regressors*b + resid_b; % Generating new data with resid_b
  b_b = regress(y_b, regressors);
  b_b_matrix(:, i) = b_b;
end

% Creating IRFs using the estimators from the bootstrapping process
% Note that the variable irf will be overwritten, but that's okay. We
% just care about cirf_b. The original cirf will not be changed.
cirf_b_matrix = zeros(h + 1, B);
for k = 1:B
  b_1_b = b_b_matrix(2:end, k);
  
  % For a moving average representation, the estimators are the IRFs.
  irf_b = b_1_b;
  
  cirf_b = zeros(h, 1);
  for i = 1:h
    for j = 1:i
      cirf_b(i) = cirf_b(i) + irf_b(j);
    end
  end
  
  a = zeros(1, 1);
  cirf_b = vertcat(a, cirf_b);
  clear a;
  
  cirf_b_matrix(:, k) = cirf_b;
end

% Gathering IRF bootstrapped SEs into a vector
cirf_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_upper = cirf + 1.96*cirf_se;
cirf_lower = cirf - 1.96*cirf_se;
%==========================================================================

%==========================================================================
%% Creating and saving cumulative IRF plot
x = linspace(0, h, length(cirf));
plot(x,cirf, 'Color', [0,0,0]);
hold on
plot(x, cirf_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Months after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Output');
legend('Impulse Response', '95% CI (Upper)', '95% CI (Lower)', 'Location', 'Southwest');

saveas(gcf, 'path\to\graphics\1b_plot.png');
close(gcf);
%==========================================================================