%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% ECO384G Problem Set 1, 1e
% Paul Le Tran, plt377
% 1 September, 2021
%==========================================================================

%==========================================================================
% Note: IRFs were created manually using 1, 18, and 36 lags of shocks. The
% number of lags can be changed in the m variable found below.
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
%ip_sa = data1(564:936, 1); % Included 1965m12 to have 1966m1 in diffs.
%chg_ln_ip_sa = diff(log(ip_sa));
% THIS IS IF YOU ARE USING ROMER^2 (2004)'S DATA INSTEAD.
chg_ln_ip_sa = data2(1:end, 5);
mp_shock = data2(1:end, 1);

% Deleting datasets as no longer needed.
clear data1 data2;
%==========================================================================

%==========================================================================
%% Making lags for each series
cd(home_dir);

% Because the purpose of this problem is to vary the number of lags for
% MP shock, we will make the lag chosen to be coded more dynamically.
% m = number of lags for MP shock
m = 1;

% Matrices include current and lagged series
% This uses's Coibion's makelags.m
chg_ln_ip_sal = makelags(chg_ln_ip_sa, 24);
mp_shockl = makelags(mp_shock, m);
%==========================================================================

%==========================================================================
%% Making every series follow the same sample period: 1970m1 - 1996m12.
% THIS IS IF YOU ARE USING FRED'S IP DATA.
%chg_ln_ip_sal = chg_ln_ip_sal(25:end, :);
% THIS IS IF YOU ARE USING ROMER^2 (2004)'S DATA INSTEAD.
chg_ln_ip_sal = chg_ln_ip_sal(25:end, :);
mp_shockl = mp_shockl((48 - m + 1):end, :);

% Making vector of current chg_ln_ip_sa for regressand
chg_ln_ip_sa = chg_ln_ip_sal(:, 1);

% Removing current vectors from lag matrices
chg_ln_ip_sal = chg_ln_ip_sal(:, 2:end);
mp_shockl = mp_shockl(:, 2:end);
%==========================================================================

%==========================================================================
%% Creating dummy variables for months (making them 324 months long)
mon = zeros(324,1);
mon(:,1) = mod(0:323, 12)+1;
dummy_data = dummyvar(mon);

% Still learning loops in matlab because I'm dumb lol.
d1 = dummy_data(:, 1);
d2 = dummy_data(:, 2);
d3 = dummy_data(:, 3);
d4 = dummy_data(:, 4);
d5 = dummy_data(:, 5);
d6 = dummy_data(:, 6);
d7 = dummy_data(:, 7);
d8 = dummy_data(:, 8);
d9 = dummy_data(:, 9);
d10 = dummy_data(:, 10);
d11 = dummy_data(:, 11);
d12 = dummy_data(:, 12);
clear mon dummy_data;

% Putting monthly dummy vars into a matrix
% Not including d12 because regression constant supposed to capture it.
dumvars = [d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11];
clear d1 d2 d3 d4 d5 d6 d7 d8 d9 10 d11 d12;
%==========================================================================

%==========================================================================
%% Creating simple regression following Romer^2 (2004) (2)
regressors = [ones(size(chg_ln_ip_sa)) dumvars chg_ln_ip_sal mp_shockl];
b = regress(chg_ln_ip_sa, regressors);

% b1 = constants + coefficients of dumvars
b_1 = b(1:12, :);

% b_2 = coefficients of chg_ln_ip_sal
b_2 = b(13:36, :);

% b_3 = coefficients of mp_shockl
b_3 = b(37:end, :);
%==========================================================================

%==========================================================================
%% Creating IRF.
% To replicate Romer^2 (2004)'s results, the mp_shock is set to 1. This
% makes the responses only the coefficients and chg_ln_ip_sa.
% Recall that we only have m lags of mp_shock, 24 lags of chg_ln_ip_sa.
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= m % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 24
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % chg_ln_ip_sa contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 25+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 24
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % chg_ln_ip_sa contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
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
b_b_matrix = zeros(length(b), B);
for i = 1:B
  resid_b = datasample(resid, length(resid)); % Reshuffling residuals
  y_b = regressors*b + resid_b; % Generating new data with resid_b
  b_b = regress(y_b, regressors);
  b_b_matrix(:, i) = b_b;
end


% Creating IRFs using the estimators from the bootstrapping process
% Note that the variable irf will be overwritten, but that's okay. We
% just care about cirf_b. The original cirf will not be changed.
cirf_b_matrix = zeros(h+1, B);
for k = 1:B
  b_2_b = b_b_matrix(13:36, k);
  b_3_b = b_b_matrix(37:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= m % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 24
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % chg_ln_ip_sa contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 25+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 24
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % chg_ln_ip_sa contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    end
  end

  % Calculating cumulative IRs for each horizon.
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

% Gathering IRF bootrstrapped SEs into a vector
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

saveas(gcf, 'path\to\graphics\1e_plot_1lags.png');
close(gcf);
%==========================================================================