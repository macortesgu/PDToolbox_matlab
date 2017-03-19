% Test for mechanism with social incentives (opinion dynamics), with Bullo stubbornnes model. 
% Replicator dynamics used for resource allocation optimization.

clear

clc

close all
tic
global G beta_ef N T_ hybrid b op devices

devices = xls2struct('devices.xlsx'); %cargar los aparatos que tienen las casas de cada usuario.
%load('Prefs.mat');
%number of populations (consumers, households)
P = 4;
% number of agents per population (devices)
N = length(devices.id);

%step_1_graph % graph network generation, for the user's opinions interaction

%step_2_opdyn % opinion dynamics



%% Users preferences
T_ = 7; %weekly period with hourly resolution.
Tperf = 10080;
 PrefsOrig = zeros(N,Tperf,P);%randi([0 1], N,T_,P);
 Perfiles = Perfcas2(P,1);

for i=1:length(devices.id)
    %update default power definition for each device with data in the
    %profiles
    devices.power(i) = ceil(max(max(Perfiles.(devices.id{i}))));
    
    %copy Perfiles to variable Prefs, changing organisation
    
    for u = 1:P
        PrefsOrig(i,:,u) = Perfiles.(devices.id{i})(:,u);
    end
    
end

%Undersampling of the profiles, to match game runtime periods.
factor = Tperf/T_;
Prefs = pUndersample(PrefsOrig, factor);
%save('Prefs','Prefs','PrefsOrig','Perfiles','devices');

%length(nonzeros(PrefsOrig(1,:,1)))
%length(nonzeros(Prefs(1,:,1)))*factor


%% Stuff

% parametros de la func. de costo agregado
%beta_ef = 450; 
b = 0;

% number of pure strategies per population
n = 2; %T_+1;%(24*7)+1;%25


%% discrete

% mass of the population
%m = 1;

% initial condition
x0 = 0; %[0.2 0.7 0.1]; 

% simulation parameters
iterations = 100;

% structure with the parameters of the game
G = struct('P', P, 'N', N, 'n', n, 'pop_wise', 0, 'f', @fitness_user_finite, 'time', iterations, 'eta', 0.02, 'revision_protocol', @pairwise_comparison);

G.R = 1;% There are as many revisions per iteration as agents are in a population.


% verify data of the game
G = definition(G);

G.T_ = T_; %assign run periods to the ciclic game.
G.Prefs = Prefs; %users' preferences for each period.

Qm = zeros(T_,1);


for run = 1:G.T_
G.period = run; %present period
G.run_finite();

Qm(run) = sum(G.Q(end,:))/1000;
end
%G.graph()
%G.graph_evolution()
%G.graph_fitness()

%% continues standard


Q_e = reshape(G.Q(end,:), [G.N,G.P]);
Q_s = reshape(G.Q(1,:), [G.N,G.P]);

user = 4;
FitEvo = zeros(G.time,G.N);
for t = 1:G.time
FitEvo(t,1:G.N) = G.F(t,user*(1:G.N));
end

figure(10)
tfit = 1:G.time;
plot(tfit,FitEvo)


% extract matrix of strategies. 
%MC. A strategy to make a consumption in a specific time t.
%x_n = vec2mat(G.X(end, :), n);

labelperiods  = periodClass(T_);

figure(1)
timem = 1:T_;
plot(timem,Qm)
xlabel(labelperiods);
ylabel('Total consumption in the society [kWh]')

%figure(2)


figure(4); 
clf
subplot(2,1,1), bar(Q_e,'DisplayName','Q_e')
titPower = strcat('Final Power allocation per user per device [W]');%, % addcost: ', num2str(addcost));
%title(titPower)
xlabel('Device') % x-axis label
ylabel(titPower) % y-axis label
legend('show', 'Location','northeastoutside')

%figure(5); 
%clf
subplot(2,1,2), bar(Q_s,'DisplayName','Q_s')
titPower = strcat('Starting Power allocation per user per device [W]');%, % addcost: ', num2str(addcost));
%title(titPower)
xlabel('Device') % x-axis label
ylabel(titPower) % y-axis label
legend('show', 'Location','northeastoutside')

%nameP = strcat('power',netType,num2str(P),social_incentives);
%savefig(4, nameP, 'compact');


%disp(['Energy use ratio: ' num2str(total_resources_alloc/sum(m))]);
timeQ = 1:size(G.Q,1);
Qavg = zeros(1,size(G.Q,1));
for i= timeQ
     Qavg(i) = sum(G.Q(i,:))/1000;
end

figure(6)
plot(timeQ,Qavg)
xlabel('iterations');
ylabel('Total consumption in the society [kWh]')


% figure(7)
% iter= 1:P;
% plot(iter,theta, iter, w, iter, susceptibility)
% legend('susceptibility','self-conf','i.c.','Location','NorthEastOutside')
% xlabel('Users');
% ylabel('Parameter value')
% savefig(7, 'parameters', 'compact');
disp('Coincide la potencia optima(prefs) con la consumida:')
PowOpt = zeros(G.P,T_);
for u = 1:G.P
    PowOpt(u,:)=squeeze(sum(power_state(Prefs(:,:,u)+1)));
end
eq(sum(PowOpt)', Qm*1000)

PotConsumida = Qm*1000
PotOptima = sum(PowOpt)'

toc
