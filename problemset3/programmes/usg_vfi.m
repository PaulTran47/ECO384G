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
% this programme computes via value function iteration (VTI) the
% equilibrium of the small open economy real business cycle (SOE RBC) model
% studied in chapter ``The Open Economy Real Business Cycle Model" of
% section ``Inducing Stationarity Through Impaticence and Global Solutions"
% in Uribe's and Schmitt-Grohe's book ``Open Economy Macroeconomic".
%=========
% END NOTE
%=========

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';

cd(home_dir);
%==========================================================================

%==========================================================================
%% Value function iteration
% Setting up calibration matrix whose columns are:
% beta dlower dupper klower kupper
% Row 1 represents the natural-debt-limit calibration
% Row 2 represents the baseline calibration
% Row 3 represents the high patience calibration
calibration =[0.9540 7.45   9.95   2.85   3.7461
              0.9540 0.6667 1.0000 2.24   3.6300
              0.9600 0      1.0000 2.5000 4.0000
]

% Performing VFI for all calibration matrix rows
for j = 1:3
  ncal = j; %pick a row of calibration (baseline 2)
  beta = calibration(ncal,1);
  dlower = calibration(ncal,2);
  dupper = calibration(ncal,3);
  klower = calibration(ncal,4);
  kupper = calibration(ncal,5);

  format compact
  load usg_tpm.mat zgrid pai

  % Setting number of grid points for natural log of technology shock
  nz = numel(zgrid);
  for i=1:nz
    a = pai(i,:);
    [~,k] = max(a);
    a(k) = 0;
    pai(i,k) = 1-sum(a);
  end

  % Calibrating more variables (as specified by Andres's problem set)
  % Degree of nominal wage rigidity
  omega = 1.455;
  
  % SS debt
  dbar = 1;
  
  % Capital adjustment cost parameter
  phi = 0.028;
  
  % Setting number of grid points for debt, d_t
  nd = 70;
  
  % Setting number of grid points for capital
  nk = 30;
  
  % Capital semielasticity of output
  alpha = 0.32;
  
  % Depreciation rate
  delta = 0.1;
  
  % Elasticity of intertemporal substitution
  sigg = 2;
  
  % World interest rate
  r = 0.04;
  
  % Discount factor
  % beta = 0.954;
  
  % Constant government spending
  g = 0;
  
  %Setting grid points for natural borrowing limit
  %dupper = 2;
  %dlower = 1;
  dgrid =  linspace(dlower,dupper,nd)';

  % SS capital stock
  %=====
  % NOTE
  %=====
  % This is obtained by running Uribe's programme edeir_ss.m
  %=========
  % END NOTE
  %=========
  kss = 3.397685279738449;
  
  % Setting grid points for capital stock
  %kupper = 4;
  %klower = 2.5
  kgrid = linspace(klower,kupper,nk)';

  % Creating d and k grids for trials
  kkgrid = repmat(kgrid',nd,1);
  kkgrid = kkgrid(:);
  ddgrid = repmat(dgrid,nk,1);

  n = nz*nd*nk;

  z = repmat(zgrid,nd*nk,1); 

  d = repmat(dgrid',nz,nk);
  d = d(:);

  k = repmat(kgrid',nz*nd,1);
  k = k(:);

  h  = ((1-alpha) * exp(z) .* k.^alpha).^(1/(alpha+omega-1));

  y = exp(z) .* k.^alpha .* h.^(1-alpha); 

  ctry = bsxfun(@(x,y) x+y,y-(1+r)*d-g,ddgrid');% consumption 

  ctry = ctry-bsxfun(@(x,y) y-(1-delta)*x+phi/2*(y-x).^2,k,kkgrid');

  ctry = bsxfun(@(x,y) x-y.^omega/omega,ctry,h);

  utry = -inf(n,nd*nk);
  utry(ctry>0) = (ctry(ctry>0).^(1-sigg)-1)/(1-sigg);

  clear ctry

  v = zeros(nz,nd,nk);
  vnew = v; 

  dist = 1;
  while dist>1e-8
    Ev = pai*v(:,:);
    Ev = repmat(Ev,nd*nk,1);

    [vnew(:),b] = max(utry+beta*Ev,[],2);

    dist = max(abs(v(:)-vnew(:)))

    v = vnew;
  end

  kp = kkgrid(b);
  dp = ddgrid(b);

  iv = kp - (1-delta)*k + phi/2 * (kp-k).^2;

  c = y -iv -g -d*(1+r) + dp;
  tb = y-c-iv - g;
  ca = -(dp-d);
  tby = tb./y;
  cay = ca./y;

  zix = repmat((1:nz)',nd*nk,1);
  [dpix,kpix] = ind2sub([nd nk],b);

  PAI = sparse(n,n);

  for i=1:n
    m1 = sub2ind([nz nd nk],1,dpix(i),kpix(i));
    PAI(i,m1:m1+nz-1) = pai(zix(i),:);
  end
  PAI = sparse(PAI);

  uPAI = usg_uncond_distrib(PAI,1e-5);

  x = zeros(nz,nd,nk);
  x(:) = uPAI;
  paiz = sum(x(:,:),2);
  paid = sum(x,3);
  paid = sum(paid)';
  paik = sum(x,1);
  paik = sum(paik,2);
  paik=paik(:);

  clear *try
  % Saving VFI results into a .mat file
  eval(['save -v7.3  usg_vfi' num2str(ncal)  '.mat'])
end
%==========================================================================