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
% this programme computes the deterministic steady-state of the EDEIR open
% economy model presented in Chapter 4 of ``Open Economy Macroeconomics" by
% Martin Uribe and Stephanie Schmitt-Grohe.
% Schmitt-Grohe.
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
%% Setting up calibration
% Time unit is a year
SIGG = 2; %mENDOZA
DELTA = 0.1; %depreciation rate
RSTAR = 0.04; %long-run interest rate
ALPHA = 0.32; %F(k,h) = k^ALPHA h^(1-ALPHA)
OMEGA = 1.455; %Frisch ela st. from Mendoza 1991
DBAR =  0.74421765717098; %debt
PSSI = 0.11135/150; %debt elasticity of interest rate
PHI = 0.028; %capital adjustment cost
RHO = 0.42; %persistence of TFP shock
ETATILDE = 0.0129; %standard deviation of innovation to TFP shock

% Implied parameters and steady state of endogenous variables
BETA = 1/(1+RSTAR); %subjective discount factor

r = RSTAR; %interest rate

d = DBAR; %debt

KAPPA = ((1/BETA - (1-DELTA)) / ALPHA)^(1/(ALPHA-1)); %k/h

h = ((1-ALPHA)*KAPPA^ALPHA)^(1/(OMEGA -1)); 

k = KAPPA * h; %capital

kfu = k; 

output = KAPPA^ALPHA * h; %output

c = output-DELTA*k-RSTAR*DBAR;

ivv = DELTA * k; %investment

tb = output - ivv - c; %trade balance

tby = tb/output;

ca = -r*d+tb;

cay = ca/output;

a = 1; %technological factor

la = ((c - h^OMEGA/OMEGA))^(-SIGG); %marginal utility of wealth
%==========================================================================