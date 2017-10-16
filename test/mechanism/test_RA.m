% example of a game with one population. three strategies per population. and combined dynamics.

% TODO:
% review how are the transitions between dynamics
% review why rd and logit are so similar in aggregated evolution. but
% differents in incentives

% population games tool box
clear
close all

global G beta_ef alpha_ef time_on time_off N T_ hybrid b q_min valCoeff addcost op



N = 100;

addcost = 0;%percent

%% opinion  dynamics parameters

%N = 600;
K = 5;
beta = 0.6; %link probability
tMax = 10e3;
c_eps = 1e-6;

G = WattsStrogatz(N,K,beta); %%graph with small world properties
A = adjacency(G);

s = rand(N,1);% %zeros(1,N)initial beliefs of users
opEps = 0.4;

fprintf(1,'å: %3.2f\n', opEps);   
[~, opinion] = hkLocal(A, s, opEps, tMax, c_eps, 'plot');

op = opinion(:,end);%zeros(1,T_);%

%% Everything else

% number of populations
P =N;

% mass of the populations -- allocated resources?
mp = 5*7; %5
m = ones(P, 1) * mp;

mp2 = 5*7; %6.88
m(P/2+1:P) = mp2;

%Definition of the electricity variables

% Dt = 0.2*[1   1   1    1    1    1 ...
%     1   1   1    1    1    1 ...
%     1   1   1    1    1    1 ...
%     1   1   1    1    1    1 ];

%Example consumption daily profile, constructed starting from total energy consumption
%historical data, in kWh

Dt = [0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18 ...
 0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18 ...
 0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18 ...
 0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18 ...
 0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18 ...
 0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18 ...
 0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18];

%valCoeffBase = 0.5*N*max(Dt)/(min(Dt)*9*sqrt(var(Dt)));%*2.5*mean(Dt)); % valor empirico asociado a la valoracion que dan los usuarios al recurso en una hora determinada.
valCoeff = 50;%valCoeffBase;%-2;

%Minimum energy consumption for all users during a time lapse.
q_min=0;%.12

pt = Dt./max(Dt)*valCoeff; %TODO: Analize pt and effect of coefficients in population stable state.

% number of strategies
T_ = length(Dt); % should be any value and work normal as well.

% valuation parameters of all agents
% valoracion heterogenea entre poblaciones
alpha_ef = zeros(N,T_);
for i=1:N/2
    alpha_ef(i,:) = pt(1:T_);%*(1+.5*i/N*0) + 0.*rand(1.T_);
end

pt = Dt./max(Dt)*(valCoeff+2);

for i=(N/2+1):N
    alpha_ef(i,:) = pt(1:T_);%*(1+.5*i/N*0) + 0.*rand(1.T_);
end

% parametros de la func. de costo agregado
beta_ef = 1;
b = 0;

% Time of the activation of either incentives or attacks IS NOT doing
% anything in the code
%time_on = 2;
%time_off = 4;


% number of pure strategies per population
n = T_+1;%(24*7)+1;%25 % MC. Energia . o potencia? que se distribuye entre los agentes. por cada poblacion.

% simulation parameters
time = 20;

% initial condition
pot = ones(N,T_+1)/(T_+1);
x0 = pot;

% structure with the parameters of the game
G = struct('P', P, 'n', n, 'f', @fitness_user, 'ode', 'ode45', 'time', time, 'tol', 0.00001, 'x0', x0, 'm', m);

% random initial condition
%G = struct('P', P, 'n', n, 'f', @fitness_user, 'ode', 'ode45', 'time', time, 'tol', 0.000001, 'm', m);

% verify data of the game
G = definition(G);

G.step = .01;
G.eta = .02;

% run dynamics
G.dynamics = {'rd'};
G.run()
T_rd = G.T;
X_rd = G.X;

% G.dynamics = {'bnn'};
% G.run()
% T_bnn = G.T;
% X_bnn = G.X;
% 
% 
% G.dynamics = {'smith'};
% G.run()
% T_smith = G.T;
% X_smith = G.X;


% % extract matrix of strategies
% %n = max(G.S);
% x_n = vec2mat(X_dyn(end. :). n);
% x = zeros(G.P. n);
% 
% for p = 1 : G.P
%     x(p. :) = x_n(p. :) * G.m(p);
% end
% 
% U = utility(x);
% 
%pause

% G.dynamics = {'logit'};
% G.run()
% T_logit = G.T;
% X_logit = G.X;



% extract matrix of strategies. 
%MC. A strategy to make a consumption in a specific time t.
x_n = vec2mat(G.X(end, :), n);
x = zeros(G.P, n);
for p = 1 : G.P
    x(p, :) = x_n(p, :) * G.m(p);
end
U = utility(x);

figure(3);
clf
plot(1:1:T_, U(1,1:T_))
tit = strcat('Utility per population, % addcost: ', num2str(addcost));
title(tit)
name = strcat('utility','%',num2str(addcost));
savefig(3, name, 'compact');

figure(4); 
clf
plot(1:1:T_, x(1,1:T_))
titPower = strcat('Power allocation per user [kW]');%, % addcost: ', num2str(addcost));
title(titPower)
xlabel('Time [h]') % x-axis label
%ylabel('sine and cosine values') % y-axis label
nameP = strcat('power','%',num2str(addcost));
hold on
plot(1:1:T_, x(10,1:T_))
savefig(4, nameP, 'compact');


%%Resource allocation verification test
for count = 1 : G.P
    resources_alloc_pu(count,1) = sum(x_n(count,1:T_));
    resources_allocated(count,1)= sum(x(count,1:T_));
end

resources_alloc_pu
resources_allocated


%if min(X_logit == )

%min(X_logit)

%graph_incentives_evolution

%G.graph_evolution()

%G.graph_state()

%sum(G.x0(1:24, :))
