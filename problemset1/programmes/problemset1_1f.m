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

% ECO384G Problem Set 1, 1f
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
data1 = xlsread(append(data_dir, '\RomerandRomerDataAppendix.xls'), 'DATA BY MONTH');
% The following data sources all have the sample period of 1965q4 -
% 1996q4. The exception is wage, which is 1979q1 - 1996q4. Furthermore,
% residential and non-residential investment are nominal because the BEA
% real series only begin in 2002.

% Real GDP
data2 = xlsread(append(data_dir, '\GDPC1.xls'));
% Real total real consumption
data3 = xlsread(append(data_dir, '\PCECC96.xls'));
% Real Real durables consumption
data4 = xlsread(append(data_dir, '\PCDG.xls'));
% Real total investment
data5 = xlsread(append(data_dir, '\GPDIC1.xls'));
% Real residential investment
data6 = xlsread(append(data_dir, '\PRFI.xls'));
% Real non-residential investment
data7 = xlsread(append(data_dir, '\PNFI.xls'));
% Real wages
data8 = xlsread(append(data_dir, '\WAGE.xls'));
% CPI
data9 = xlsread(append(data_dir, '\CPI.xls'));
% PPI
data10 = xlsread(append(data_dir, '\PPI.xls'));
%==========================================================================

%==========================================================================
%% Creating time series
% Converting MP shock series into quarterly frequency
mp_shock_m = data1(1:end, 1);
% 1966q1 - 1996q4 = 124 qtrs
mp_shock = zeros(124, 1);
for i = 1:length(mp_shock)
  % Creating index that represents starting month to calculate qtr value
  start_m = 3*i - 2;
  start_m;
  avg = (mp_shock_m(start_m) + mp_shock_m(start_m + 1) + mp_shock_m(start_m + 1))/3
  mp_shock(i) = avg;
end

% Creating natural logs and difference of natural logs. Latter will all be
% 1966q1 - 1996q4, except for wages.
gdp = data2(1:end, 1);
ln_gdp = log(gdp);
chg_ln_gdp = diff(ln_gdp);

pce = data3(1:end, 1);
ln_pce = log(pce);
chg_ln_pce = diff(ln_pce);

pced = data4(1:end, 1);
ln_pced = log(pced);
chg_ln_pced = diff(ln_pced);

bfi = data5(1:end, 1);
ln_bfi = log(bfi);
chg_ln_bfi = diff(ln_bfi);

bfir = data6(1:end, 1);
ln_bfir = log(bfir);
chg_ln_bfir = diff(ln_bfir);

bfin = data7(1:end, 1);
ln_bfin = log(bfin);
chg_ln_bfin = diff(ln_bfin);

% The diff in log(wages) will begin in 1979q2.
wage = data8(1:end, 1);
ln_wage = log(wage);
chg_ln_wage = diff(ln_wage);

cpi = data9(1:end, 1);
ln_cpi = log(cpi);
chg_ln_cpi = diff(ln_cpi);

ppi = data10(1:end, 1);
ln_ppi = log(ppi);
chg_ln_ppi = diff(ln_ppi);

% Deleting datasets as no longer needed.
clear data1 data2 data3 data4 data5 data6 data7 data8 data9 data10;
%==========================================================================

%==========================================================================
%% Making lags for each series
cd(home_dir);

% Using 12 lags for MP shock
mp_shockl = makelags(mp_shock, 12);
% Making lagged matrix follow sample period of 1970q1 - 1996q4.
mp_shockl = mp_shockl(5:end, :);

% Using 8 lags for every other series
% Making lagged matrix follow sample period of 1970q1 - 1996q4, except for
% wages.
chg_ln_gdpl = makelags(chg_ln_gdp, 8);
chg_ln_gdpl = chg_ln_gdpl(9:end, :);

chg_ln_pcel = makelags(chg_ln_pce, 8);
chg_ln_pcel = chg_ln_pcel(9:end, :);

chg_ln_pcedl = makelags(chg_ln_pced, 8);
chg_ln_pcedl = chg_ln_pcedl(9:end, :);

chg_ln_bfil = makelags(chg_ln_bfi, 8);
chg_ln_bfil = chg_ln_bfil(9:end, :);

chg_ln_bfirl = makelags(chg_ln_bfir, 8);
chg_ln_bfirl = chg_ln_bfirl(9:end, :);

chg_ln_bfinl = makelags(chg_ln_bfin, 8);
chg_ln_bfinl = chg_ln_bfinl(9:end, :);

chg_ln_cpil = makelags(chg_ln_cpi, 8);
chg_ln_cpil = chg_ln_cpil(9:end, :);

chg_ln_ppil = makelags(chg_ln_ppi, 8);
chg_ln_ppil = chg_ln_ppil(9:end, :);

% Making regressand vectors
chg_ln_gdp = chg_ln_gdpl(:, 1);
chg_ln_pce = chg_ln_pcel(:, 1);
chg_ln_pced = chg_ln_pcedl(:, 1);
chg_ln_bfi = chg_ln_bfil(:, 1);
chg_ln_bfir = chg_ln_bfirl(:, 1);
chg_ln_bfin = chg_ln_bfinl(:, 1);
chg_ln_cpi = chg_ln_cpil(:, 1);
chg_ln_ppi = chg_ln_ppil(:, 1);

% Removing current vectors from lag matrices
mp_shockl = mp_shockl(:, 2:end);
chg_ln_gdpl = chg_ln_gdpl(:, 2:end);
chg_ln_pcel = chg_ln_pcel(:, 2:end);
chg_ln_pcedl = chg_ln_pcedl(:, 2:end);
chg_ln_bfil = chg_ln_bfil(:, 2:end);
chg_ln_bfirl = chg_ln_bfirl(:, 2:end);
chg_ln_bfinl = chg_ln_bfinl(:, 2:end);
chg_ln_cpil = chg_ln_cpil(:, 2:end);
chg_ln_ppil = chg_ln_ppil(:, 2:end);
%==========================================================================

%==========================================================================
%% Running OLS regressions and creating IRFs for each regressand.
%========
% y = GDP
%========
% Creating regressor matrix
regressors = [ones(size(chg_ln_gdp)) chg_ln_gdpl mp_shockl];

% Running OLS regression
b = regress(chg_ln_gdp, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_gdp = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_gdp(i) = cirf_gdp(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_gdp = vertcat(a, cirf_gdp);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_gdp - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_gdp_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_gdp_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_gdp_upper = cirf_gdp + 1.96*cirf_gdp_se;
cirf_gdp_lower = cirf_gdp - 1.96*cirf_gdp_se;

%========
% y = PCE
%========
% Creating regressor matrix
regressors = [ones(size(chg_ln_pce)) chg_ln_pcel mp_shockl];

% Running OLS regression
b = regress(chg_ln_pce, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_pce = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_pce(i) = cirf_pce(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_pce = vertcat(a, cirf_pce);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_pce - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_pce_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_pce_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_pce_upper = cirf_pce + 1.96*cirf_pce_se;
cirf_pce_lower = cirf_pce - 1.96*cirf_pce_se;

%=========
% y = PCED
%=========
% Creating regressor matrix
regressors = [ones(size(chg_ln_pced)) chg_ln_pcedl mp_shockl];

% Running OLS regression
b = regress(chg_ln_pced, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_pced = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_pced(i) = cirf_pced(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_pced = vertcat(a, cirf_pced);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_pced - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_pced_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_pced_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_pced_upper = cirf_pced + 1.96*cirf_pced_se;
cirf_pced_lower = cirf_pced - 1.96*cirf_pced_se;

%========
% y = BFI
%========
% Creating regressor matrix
regressors = [ones(size(chg_ln_bfi)) chg_ln_bfil mp_shockl];

% Running OLS regression
b = regress(chg_ln_bfi, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_bfi = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_bfi(i) = cirf_bfi(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_bfi = vertcat(a, cirf_bfi);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_bfi - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_bfi_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_bfi_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_bfi_upper = cirf_bfi + 1.96*cirf_bfi_se;
cirf_bfi_lower = cirf_bfi - 1.96*cirf_bfi_se;

%=========
% y = BFIR
%=========
% Creating regressor matrix
regressors = [ones(size(chg_ln_bfir)) chg_ln_bfirl mp_shockl];

% Running OLS regression
b = regress(chg_ln_bfir, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_bfir = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_bfir(i) = cirf_bfir(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_bfir = vertcat(a, cirf_bfir);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_bfir - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_bfir_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_bfir_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_bfir_upper = cirf_bfir + 1.96*cirf_bfir_se;
cirf_bfir_lower = cirf_bfir - 1.96*cirf_bfir_se;

%=========
% y = BFIN
%=========
% Creating regressor matrix
regressors = [ones(size(chg_ln_bfin)) chg_ln_bfinl mp_shockl];

% Running OLS regression
b = regress(chg_ln_bfin, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_bfin = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_bfin(i) = cirf_bfin(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_bfin = vertcat(a, cirf_bfin);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_bfin - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_bfin_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_bfin_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_bfin_upper = cirf_bfin + 1.96*cirf_bfin_se;
cirf_bfin_lower = cirf_bfin - 1.96*cirf_bfin_se;

%========
% y = CPI
%========
% Creating regressor matrix
regressors = [ones(size(chg_ln_cpi)) chg_ln_cpil mp_shockl];

% Running OLS regression
b = regress(chg_ln_cpi, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_cpi = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_cpi(i) = cirf_cpi(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_cpi = vertcat(a, cirf_cpi);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_cpi - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_cpi_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_cpi_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_cpi_upper = cirf_cpi + 1.96*cirf_cpi_se;
cirf_cpi_lower = cirf_cpi - 1.96*cirf_cpi_se;

%========
% y = PPI
%========
% Creating regressor matrix
regressors = [ones(size(chg_ln_ppi)) chg_ln_ppil mp_shockl];

% Running OLS regression
b = regress(chg_ln_ppi, regressors);

% Constant
b_1 = b(1, 1);

% Coefficients associated with lagged regressand matrix
b_2 = b(2:9, 1);

% Coefficients associated with MP shock
b_3 = b(10:end, 1);

% Creating IRF
h = 48
irf = zeros(h, 1);
irf(1) = b_3(1);

% Calculating individual IRs for each horizon.
for i = 2:h
  if i <= 12 % Refers to number of mp_shock lags
    irf(i) = b_3(i); % mp_shock contribution
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  else % No more contribution from mp_shock for horizons 13+
    irf(i) = 0;
    for j = 1:(i - 1)
      if i - j <= 8
        irf(i) = irf(i) + b_2(i - j)*irf(j);  % regressand contributions
      else
        irf(i) = irf(i);
      end
    end
  end
end

% Calculating cumulative IRs for each horizon.
cirf_ppi = zeros(h, 1);
for i = 1:h
  for j = 1:i
    cirf_ppi(i) = cirf_ppi(i) + irf(j);
  end
end

a = zeros(1, 1);
cirf_ppi = vertcat(a, cirf_ppi);
clear a;

% Obtaining bootstrapped SEs
% Will do the bootstrapping process B times
B = 10000;

% Obtain the original regression residuals
y_hat = regressors*b;
resid = chg_ln_ppi - y_hat;

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
  b_2_b = b_b_matrix(2:9, k);
  b_3_b = b_b_matrix(10:end, k);
    
  irf_b = zeros(h, 1);
  irf_b(1) = b_3_b(1);
  
  % Calculating individual IRs for each horizon.
  for i = 2:h
    if i <= 12 % Refers to number of mp_shock lags
      irf_b(i) = b_3_b(i); % mp_shock contribution
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
        else
          irf_b(i) = irf_b(i);
        end
      end
    else % No more contribution from mp_shock for horizons 13+
      irf_b(i) = 0;
      for j = 1:(i - 1)
        if i - j <= 8
          irf_b(i) = irf_b(i) + b_2(i - j)*irf_b(j);  % regressand contributions
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
cirf_ppi_se = zeros(h + 1, 1);
for i = 1:(h + 1)
  cirf_ppi_se(i, 1) = std(cirf_b_matrix(i, :)');
end

% Constructing 95% CIs
cirf_ppi_upper = cirf_ppi + 1.96*cirf_ppi_se;
cirf_ppi_lower = cirf_ppi - 1.96*cirf_ppi_se;
%==========================================================================

%==========================================================================
%% Creating and saving cumulative IRF plot
x = linspace(0, h, length(cirf_gdp));
tiledlayout(4, 2);

% GDP
nexttile
plot(x,cirf_gdp, 'Color', [0,0,0]);
hold on
plot(x, cirf_gdp_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_gdp_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Real GDP');

% PCE
nexttile
plot(x,cirf_pce, 'Color', [0,0,0]);
hold on
plot(x, cirf_pce_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_pce_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Real PCE');

% PCED
nexttile
plot(x,cirf_pced, 'Color', [0,0,0]);
hold on
plot(x, cirf_pced_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_pced_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Real PCE, Durables');

% BFI
nexttile
plot(x,cirf_bfi, 'Color', [0,0,0]);
hold on
plot(x, cirf_bfi_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_bfi_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Real BFI');

% BFIR
nexttile
plot(x,cirf_bfir, 'Color', [0,0,0]);
hold on
plot(x, cirf_bfir_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_bfir_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Real BFI, Residential');

% BFIN
nexttile
plot(x,cirf_bfin, 'Color', [0,0,0]);
hold on
plot(x, cirf_bfin_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_bfin_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of Real BFI, Non-Residential');

% CPI
nexttile
plot(x,cirf_cpi, 'Color', [0,0,0]);
hold on
plot(x, cirf_cpi_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_cpi_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of CPI');

% PPI
nexttile
plot(x,cirf_ppi, 'Color', [0,0,0]);
hold on
plot(x, cirf_ppi_upper, 'LineStyle', '--', 'Color', [0,0,0]);
hold on
plot(x, cirf_ppi_lower, 'LineStyle', '--', 'Color', [0,0,0]);
hold off
grid on
xlabel('Quarters after Shock');
ylabel('Percents in decimal form');
hold on
title('Impulse Response of PPI');

saveas(gcf, 'path\to\graphics\1f_plot.png');
close(gcf);
%==========================================================================