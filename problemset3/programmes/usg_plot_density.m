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
% this programme produces figures with ergodic debt distributions for the
% SOE RBC model solved with global methods (computed in usg_vfi.m).
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
%% Plotting debt ergodic distributions
% Plotting debt ergodic distribution when VFI performed with baseline
% calibration
load usg_vfi2
dgrid2 = dgrid;
paid2 = paid;
% subplot(2,2,[3 4])
plot(dgrid,paid,'-','linewidth',2)
h=title(['Debt ergodic distribution, $\bar{d}=$ ' num2str(dupper,1) ' and $\beta(1+r^*)=$ ' num2str(beta*(1+r),4)] , 'interpreter', 'LaTeX')
xlabel('Debt, d')
xlim([dlower dupper])

% % Plotting debt ergodic distribution when VFI performed with
% % natural-debt-limit calibration
% load usg_vfi1
% dgrid1 = dgrid;
% paid1 = paid;
% subplot(2,2,1)
% plot(dgrid,paid,'-','linewidth',2)
% h=title(['$\bar{d}=$ '  num2str(dupper,  '%0.3g')  ' and $\beta(1+r^*)=$ ' num2str(beta*(1+r),4)] ,'interpreter','LaTeX')
% xlabel('Debt, d')
% xlim([dlower dupper])
% 
% % Plotting debt ergodic distribution when VFI performed with high patience
% % calibration
% load usg_vfi3
% dgrid3 = dgrid;
% paid3 = paid;
% subplot(2,2,2)
% plot(dgrid,paid,'-','linewidth',2)
% h=title(['$\bar{d}=$ ' num2str(dupper,1) ' and $\beta(1+r^*)=$ ' num2str(beta*(1+r),4)] ,'interpreter','LaTeX')
% xlabel('Debt, d')
% xlim([dlower dupper])

shg

saveas(gcf, 'path\to\graphics\3_ergodic_d_dists.png');
close(gcf);
%==========================================================================