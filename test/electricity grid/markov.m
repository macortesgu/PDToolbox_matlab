TRANS = [.9 .1; .05 .95;];

EMIS = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6;...
7/12, 1/12, 1/12, 1/12, 1/12, 1/12];

[seq,states] = hmmgenerate(1000,TRANS,EMIS);

% Cómo escogió los valores de lambda para cada aparato? Posiblemente
% teniendo en cuenta una cantidad de energía consumida promedio acumulada.
% cómo escogió la proporción de energía para cada aparato? con los
% criterios de apertura de puertas y eso y la duración de cada evento--?
% los sucesos no son cadenas de markov, porque el estado futuro depende en
% algunas ocasiones de la historia. e.gr. cuando se prende y sigue prendido
% cierto tiempo.

devices = xls2struct('devices.xlsx'); %cargar los aparatos que tienen las casas de cada usuario.
%load('Prefs.mat');
%number of populations (consumers, households)
P = 1;
% number of agents per population (devices)
N = length(devices.id);

%step_1_graph(P) % graph network generation, for the user's opinions interaction
%step_2_opdyn % opinion dynamics

Users preferences
T_ = 14;%168/4;%10080/5; %weekly period with hourly resolution.
Tperf = 10080; % Perfiles generados con resolucion de 1 min.

 PrefsOrig = zeros(N,Tperf,P);%randi([0 1], N,T_,P);
 Perfiles = Perfcas2(P,1);