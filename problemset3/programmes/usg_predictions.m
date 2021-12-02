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
% this programme produces a table with unconditional second moments for the
% SOE RBC model with impatience solved with global methods (computed in
% usg_vfi.m).
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
%% Simulating model to produce table with unconditional second moments
format compact

%=====
% NOTE
%=====
% The files usg_vfi1.mat, usg_vfi2.mat, and usg_vfi3.mat used in this
% programme must be  produced by running usg_vfi.m in and setting mcal to
% the values specified as follows:
% usg_vfi1.mat contains policy functions for the Natural-Debt-Limit
% Calibration (set ncal = 1 in usg_vfi.m). 
% usg_vfi2.mat contains policy functions for the Baseline calibration
% (set ncal = 2 in usg_vfi.m). And
% usg_vfi3.mat contains policy functions for the  High-Paticence
% calibration (set ncal = 3 in usg_vfi.m). 
%=========
% END NOTE
%=========
filename = 'usg_vfi2'
eval(['load ' filename])

%unconditional means of y c iv h
x = [y c iv h tby cay];

%unconditional expectations
Ex = uPAI'*x;

%divide y c iv h by their respective uncond. means (as opposed to taking logs, because iv is negative in some states)
x1 = [[y/Ex(1) c/Ex(2) iv/Ex(3) h/Ex(4)-1] tby-Ex(5) cay-Ex(6)];

%unconditional first and second moments
[Ex1,Vx1,SDx1,Corrx1,Scorrx1] = usg_moments_tpm(PAI,x1,1,uPAI);

%extract sedond moments for table (std dev, autocorr, and corr with output)

disp('Columns are: std, serial correlation, corr. w. output.  Rows are: y, c, i, h, tb/y, ca/y')
usg_table = [SDx1'*100 Scorrx1' Corrx1(:,1)]

eval(['save predictions_' filename ' beta r  dupper Ex SDx1 '])
%==========================================================================