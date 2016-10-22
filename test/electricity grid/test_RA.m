% example of a game with one population. three strategies per population. and combined dynamics.

% TODO:
% review how are the transitions between dynamics
% review why rd and logit are so similar in aggregated evolution. but
% differents in incentives

% population games tool box
clear
close all

global G beta_ef alpha_ef time_on time_off mp N T_ hybrid b q_min Dt ptCoeff

q_min=0;

N = 2;

%Definition of the electricity variables
%TODO: Include here variation with constructed residential user' load
%curve. having the total energy limit.
%Dt = 15*[51.8743   50.0011   48.6104    48.6384    51.1276    58.7756 ...
%    61.0654   65.0167   69.6593    71.6363    75.3904    76.2807 ...
%    73.4635   73.3627   74.6492    75.1194    74.8689    74.1951 ...
%    78.2569   85.8935   83.5392    77.9073    68.6800   60.5177];



% Dt = [1   1   1    1    1    1 ...
%     1   1   1    1    1    1 ...
%     1   1   1    1    1    1 ...
%     1   1   1    1    1    1 ];

pt = Dt./max(Dt)*ptCoeff; %TODO: Analize pt and effect of coefficients in population stable state.

% number of strategies
T_ = length(Dt);

% valuation parameters of all agents
% valoracion homogenea entre poblaciones
alpha_ef = zeros(N,T_);
for i=1:N
    alpha_ef(i,:) = pt(1:T_);%*(1+.5*i/N*0) + 0.*rand(1.T_);
end

% parametros de la func. de costo agregado
beta_ef = 1;
b = 0;

% Time of the activation of either incentives or attacks
time_on = 2;
time_off = 4;


% number of populations
P =N;

% number of pure strategies per population
n = 24;%25 % MC. Energia . o potencia? que se distribuye entre los agentes. por cada poblacion.

% mass of the populations -- allocated resources?
mp = 6.8825; %14
m = ones(P, 1) * mp;

% simulation parameters
time = 70;

% initial condition
pot = ones(N,T_)/(T_);%pot = ones(N.T_+1)/(T_+1);
x0 = pot;

% structure with the parameters of the game
G = struct('P', P, 'n', n, 'f', @fitness_user, 'ode', 'ode113', 'time', time, 'tol', 0.00001, 'x0', x0, 'm', m);

% random initial condition
%G = struct('P'. P. 'n'. n. 'f'. @fitness_user. 'ode'. 'ode45'. 'time'. time. 'tol'. 0.000001. 'm'. m);

% verify data of the game
G = definition(G);

G.step = .01;
G.eta = .02;

% run different dynamics
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
%MC. A strategy should mean a specific consumption in a specific time t.
%Apparently, what I am getting now is that strategies are consumption
%intervals (i.e. hours, 24.)
x_n = vec2mat(G.X(end, :), n);
x = zeros(G.P, n);
for p = 1 : G.P
    x(p, :) = x_n(p, :) * G.m(p);
end
U = utility(x);

figure(3);
clf
plot(1:1:T_, U)
tit = strcat('Utility per population, ptCoeff: ', num2str(ptCoeff));
title(tit)
name = strcat('utility','pt',num2str(ptCoeff));
savefig(3, name, 'compact');

figure(4); 
clf
plot(1:1:T_, x)
titPower = strcat('Power allocation per population, ptCoeff: ', num2str(ptCoeff));
title(titPower)
nameP = strcat('power','pt',num2str(ptCoeff));
savefig(4, nameP, 'compact');

%%Resource allocation verification test
for count = 1 : G.P
    resources_alloc_pu(count,1) = sum(x_n(count,:));
    resources_allocated(count,1)= sum(x(count,:));
end

resources_alloc_pu
resources_allocated


%if min(X_logit == )

%min(X_logit)

%graph_incentives_evolution

%G.graph_state()

%sum(G.x0(1:24, :))
