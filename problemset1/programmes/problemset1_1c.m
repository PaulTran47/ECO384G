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

% ECO384G Problem Set 1, 1c
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
ln_ip_sa = log(ip_sa);
chg_ln_ip_sa = diff(ln_ip_sa);
mp_shock = data2(1:end, 1);

% Making MP shock follow sample period of 1970m1 - 1996m12.
mp_shock = mp_shock(49:end, :);

% Creating separate IP series for the purposes of backwards shifting
% later on for LP. The time series will go from 1970m1 - 2000m12.
ip_sa_f = data1(613:984, 1);
ln_ip_sa_f = log(ip_sa_f);

% Deleting datasets as no longer needed.
clear data1 data2;
%==========================================================================

%==========================================================================
%% Making lags for IP
cd(home_dir);
chg_ln_ip_sal = makelags(chg_ln_ip_sa, 24);
% Making matrix follow sample period of 1970m1 - 1996m12.
chg_ln_ip_sal = chg_ln_ip_sal(25:end, :);

% Making vector of current chg_ln_ip_sa for regressand
chg_ln_ip_sa = chg_ln_ip_sal(:, 1);

% Removing current vectors from lag matrix
chg_ln_ip_sal = chg_ln_ip_sal(:, 2:end);
%==========================================================================

%==========================================================================
%% Doing local projections for h = 1:48
h = 48;
irf = zeros(h, 1);
irf_se = zeros(h, 1);
% In order to get our horizon values, we will be shifting our regressand
% backwards in order to match the horizon. Essentially, for horizon h,
% what we do is make a series where y(t) is replaced with y(t + h) - y(t).

% Creating matrix where each column represents regressand y_h = y_t+h - y_t
y_h_matrix = zeros(length(chg_ln_ip_sa), h);
for i = 1:h % i represents the regressand for each horizon (column).
  y_h = zeros(length(chg_ln_ip_sa), 1);
  for j = 1:length(ln_ip_sa_f) - h % j represents the horizon itself when calculating the values in the time series
    y_h(j) = ln_ip_sa_f(j + i) - ln_ip_sa_f(j);
  end
  y_h_matrix(:, i) = y_h;
end

% Calculating individual IRs for each horizon.
irf = zeros(h, 1);
irf_se = zeros(h, 1);
irf_upper = zeros(h, 1);
irf_lower = zeros(h, 1);
for i = 1:h
  % Creating variables for OLS regression.
  regressand = y_h_matrix(:, i);
  % We are leaving vector of ones out of regressor matrix because we will
  % be using fitlm(X, y) instead of general regress(y, X).
  regressors = [chg_ln_ip_sal mp_shock];
  
  % Running LP OLS regression. Obtaining coefficients and HAC SEs.
  fit = fitlm(regressors, regressand);
  [~, se_hac, b_h] = hac(fit, 'type', 'HC', 'weights', 'HC1', 'display', 'off');
  
  % Saving MP shock coefficient (IR of output at horizon h).
  irf(i) = b_h(26, 1);
  
  % Saving MP shock coefficient HAC SE.
  irf_se(i) = se_hac(26, 1);
end

% Constructing 95% CIs.
irf_upper = irf + 1.96*irf_se;
irf_lower = irf - 1.96*irf_se;

a = zeros(1, 1);
irf = vertcat(a, irf);
irf_upper = vertcat(a, irf_upper);
irf_lower = vertcat(a, irf_lower);
% Scaling IRFs and IRF CIs by 100x so percentages are in integer form.
%irf = irf*100;
%irf_upper = irf_upper*100;
%irf_lower = irf_lower*100;
clear a;
%==========================================================================

%==========================================================================
%% Creating and saving cumulative IRF plot
x = linspace(0, 49, length(irf));
plot(x,irf, 'Color', [0,0,0]);
hold on
plot(x, irf_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, irf_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Months after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Output');
legend('Impulse Response', '95% CI (Upper)', '95% CI (Lower)', 'Location', 'Southwest');

saveas(gcf, 'path\to\graphics\1c_plot.png');
close(gcf);
%==========================================================================