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
% this programme produces impulse response functions to a technology shock
% implied by the SOE RBC model with impatience solved with global methods
% (computed in usg_vfi.m).
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
%% Producing IRFs
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
load usg_vfi2 PAI uPAI y c iv tby cay h z zix

% Initial state of the economy (zix(6) is roughly 1 SD of the innovation of
% z).
%=====
% NOTE
%=====
% s is a vector because all that is in the information set is the knowledge
% that zix = 6 in the initial period.
%=========
% END NOTE
%=========
s = find(zix == 6);

X1 = [y c iv h];
EX1 = uPAI'*X1;
x1 = bsxfun(@(x,y) (x./y-1)*100,X1,EX1);
X2 = [tby z cay];
EX2 = uPAI'*X2;
x2 = bsxfun(@(x,y) (x-y)*100,X2,EX2);
x = [x1 x2];

T = 11;
IR = usg_ir_tpm(PAI,uPAI,x,s,T);

IR = IR/IR(1,6);
t = (0:10)';

subplot(3,2,1)
plot(t,IR(:,1),'linewidth',1.5)
title('Output')
ylabel('% dev from mean')
subplot(3,2,1)

subplot(3,2,2)
plot(t,IR(:,2),'linewidth',1.5)
title('Consumption')
ylabel('% dev from mean')


subplot(3,2,3)
plot(t,IR(:,3),'linewidth',1.5)
title('Investment')
ylabel('% dev from mean')

subplot(3,2,4)
plot(t,IR(:,4),'linewidth',1.5)
title('Hours')
ylabel('% dev from mean')

subplot(3,2,5)
plot(t,IR(:,5),'linewidth',1.5)
title('Trade Balance / Output')
ylabel('dev from mean in %')


subplot(3,2,6)
plot(t,IR(:,6),'linewidth',1.5)
title('TFP Shock')
ylabel('% dev from mean')

saveas(gcf, 'path\to\graphics\3_irf.png');
close(gcf);
%==========================================================================