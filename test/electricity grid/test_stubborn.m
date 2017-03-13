% Test for mechanism with social incentives (opinion dynamics), with Bullo stubbornnes model. 
% Replicator dynamics used for resource allocation optimization.
clc
clear

close all

global G beta_ef alpha_ef time_on time_off N T_ hybrid b q_min valCoeff addcost op mp devices

%number of populations (consumers, households)
P = 10;
addcost = 0;

step_1_graph % graph network generation, for the user's opinions interaction

step_2_opdyn % opinion dynamics

devices = xls2struct('devices.xlsx'); %cargar los aparatos que tienen las casas de cada usuario.

%% Everything else
% number of agents per population (devices)
N = 6;

%Example consumption daily profile, constructed starting from total energy consumption
%historical data, in kWh
%xlRange = 'A:B';
%DtOrig = xlsread('PerfilConsumoEjemplo.xlsx','Gen',xlRange); %load weekly profile from Excel, with minute-resolution

%kWhOrig = trapz(DtOrig)/(1000*60);%original signal has minute resolution, and energy is measured in kWh
%freq_p=1;
%freq_q=5;
%d_factor = 10;

% tx = 1:1:length(DtOrig);
% tz = 1:d_factor*6:length(DtOrig);

%Dt_temp = decimate(DtOrig(:,1), d_factor, 40, 'FIR');
%Dt = decimate(Dt_temp, 6, 40, 'FIR');
%kWh = trapz(Dt)/(1000*60/(d_factor*6));

% figure();
% plot(tx,DtOrig(:,1),'+-',tz,Dt(:,1),'o:')
% legend('original','resampled')


%% Users preferences
T_ = 48; %weekly period with hourly resolution.
Prefs = randi([0 1], N,T_,G.P);

%% Stuff

% parametros de la func. de costo agregado
beta_ef = 450;
b = 0;

% number of pure strategies per population
n = 2; %T_+1;%(24*7)+1;%25


%% discrete

% mass of the population
%m = 1;

% initial condition
%x0 = [0.2 0.7 0.1]; 

% simulation parameters
iterations = 200;


% structure with the parameters of the game
G = struct('P', P, 'N', N, 'n', n, 'f', @fitness_user_finite, 'time', iterations, 'eta', 0.02, 'revision_protocol', @pairwise_comparison);

G.R = 1;

% verify data of the game
G = definition(G);



G.T_ = T_; %assign run periods to the ciclic game.
G.Prefs = Prefs; %users' preferences for each period.

Qm = zeros(T_,1);
for run = 1:G.T_
G.period = run; %first period
G.run_finite();
%G.graph()
%G.graph_evolution()
%G.graph_fitness()
Qm(run) = sum(G.Q(run,:))/1000;


end

%% continues standard


Q_e = reshape(G.Q(end,:), [G.N,G.P]);
Q_s = reshape(G.Q(1,:), [G.N,G.P]);

% verify data of the game
%G = definition(G);

%G.step = .01;
%G.eta = .02;

% run dynamics
% G.dynamics = {'rd'};
% G.run()
% T_rd = G.T;
% X_rd = G.X;



% extract matrix of strategies. 
%MC. A strategy to make a consumption in a specific time t.
%x_n = vec2mat(G.X(end, :), n);



figure(1)
timem = 1:T_;
plot(timem,Qm)
xlabel('iterations');
ylabel('Total consumption in the society [kWh]')


figure(4); 
clf
bar(Q_e,'DisplayName','Q_e')
titPower = strcat('Final Power allocation per user per device [W]');%, % addcost: ', num2str(addcost));
%title(titPower)
xlabel('Device') % x-axis label
ylabel(titPower) % y-axis label
legend('show', 'Location','northeastoutside')

figure(5); 
clf
bar(Q_s,'DisplayName','Q_s')
titPower = strcat('Starting Power allocation per user per device [W]');%, % addcost: ', num2str(addcost));
%title(titPower)
xlabel('Device') % x-axis label
ylabel(titPower) % y-axis label
legend('show', 'Location','northeastoutside')



%nameP = strcat('power',netType,num2str(P),social_incentives);
%savefig(4, nameP, 'compact');


%disp(['Energy use ratio: ' num2str(total_resources_alloc/sum(m))]);
timeQ = 1:length(G.Q);
for i= timeQ
     Qavg(i) = sum(G.Q(i,:))/1000;
end

figure(6)
plot(timeQ,Qavg)
xlabel('iterations');
ylabel('Total consumption in the society [kWh]')


figure(7)
iter= 1:P;
plot(iter,theta, iter, w, iter, susceptibility)
legend('susceptibility','self-conf','i.c.','Location','NorthEastOutside')
xlabel('Users');
ylabel('Parameter value')
savefig(7, 'parameters', 'compact');
%if min(X_logit == )

%min(X_logit)

%graph_incentives_evolution

%G.graph_evolution()

%G.graph_state()

%sum(G.x0(1:24, :))
