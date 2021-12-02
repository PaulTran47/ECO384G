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
% this programme produces a table useful for understanding the calibration
% chosen for the computation of the SOE RBC model solved with global
% methods (computed in usg_vfi.m).
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
%% Extracting several calibration parameters from tables produced in usg_predictions.m

%Baseline calibration
load predictions_usg_vfi2.mat  
%produced by running 
%usg_predictions.m in
%and setting 
%filename = 'usg_vfi2'
usg_table_targets(1,1:5) = [beta beta*(1+r) dupper Ex(5) SDx1(3)];

%Natural-debt-limit calibration
load predictions_usg_vfi1.mat  
%produced by running 
%usg_predictions.m 
% and setting 
%filename = 'usg_vfi1'
usg_table_targets(2,1:5) = [beta beta*(1+r) dupper Ex(5) SDx1(3)];

%High-patience calibration
load predictions_usg_vfi3.mat  
%produced by running 
%usg_predictions.m 
%and setting 
%filename = 'usg_vfi3'
usg_table_targets(3,1:5) = [beta beta*(1+r) dupper Ex(5) SDx1(3)];

usg_table_targets;
%==========================================================================