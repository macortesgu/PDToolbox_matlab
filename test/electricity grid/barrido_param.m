%% Barrido paramétrico

global Dt ptCoeff

Dt = [0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18];

% ptCoeff = 20;
% %pt = Dt./max(Dt)*20; %TODO: Analize pt and effect of coefficients in population stable state.
% 
% test_RA
% 
% ptCoeff = 15;
% test_RA
% 
% ptCoeff = 10;
% test_RA
% 
% ptCoeff = 30;
% test_RA

ptCoeff = 500;
test_RA

SubPlot