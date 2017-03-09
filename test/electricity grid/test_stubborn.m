% Test for mechanism with social incentives (opinion dynamics), with Bullo stubbornnes model. 
% Replicator dynamics used for resource allocation optimization.
clc
clear

close all

global G beta_ef alpha_ef time_on time_off N T_ hybrid b q_min valCoeff addcost op mp

N = 10;
addcost = 0;

step_1_graph

s = linspace(1,0.2,N);% %zeros(N,1)initial beliefs of users
tMax = 1e3; %Max. iterations for op dynamics
c_eps = 1e-4;%1e-6; error aceptable para la convergencia
%opEps = 0.4; Bounded confidence parameter for HK op dyn model

%% susceptibility (1-stubborness) theta of the agents
theta = (1-(logspace(0, 3, N)/1000))';%column vector of susceptibility
%define self confidence
w1 = logspace(0, 3, N)/1000;
w2 = fliplr(w1);
w = (w1+w2)/2;
w = w';
%fprintf(1,'å: %3.2f\n', opEps);   
%[~, opinion] = hkLocal(A, s, opEps, tMax, c_eps, 'plot');

opinion = bullo2(s,N,Aw,theta,tMax,w, c_eps, 'plot');

social_incentives = 'none';
switch social_incentives
    case 'sInc'
        op = opinion(:,end); %zeros(1,T_);%
        
    otherwise
        op = zeros(N,1); %No use of social incentives
        social_incentives = 'noneInc';
end



%save('Caso1.mat','Gw','opinion'); 

%% Everything else
% number of populations
P =N;

%Example consumption daily profile, constructed starting from total energy consumption
%historical data, in kWh
xlRange = 'A:B';
DtOrig = xlsread('PerfilConsumoEjemplo.xlsx','Gen',xlRange); %load weekly profile from Excel, with minute-resolution

kWhOrig = trapz(DtOrig)/(1000*60);%original signal has minute resolution, and energy is measured in kWh
%freq_p=1;
%freq_q=5;
%d_factor = 5;

%tx = 1:1:length(DtOrig);
%tz = 1:d_factor:length(DtOrig);

%Dt = resample(DtOrig,tx,'pchip');
Dt = DtOrig;%decimate(DtOrig(:,1), d_factor, 40, 'FIR');
kWh = kWhOrig;%trapz(Dt)/(1000*60/(freq_q/freq_p));



%figure();
%plot(tx,DtOrig(:,1),'+-',tz,Dt(:,1),'o:')
%legend('original','resampled')


% mass of the populations -- allocated resources?
mp = kWh(1);
m = ones(P, 1) * mp;

mp2 = kWh(2); 
m(P/2+1:P) = mp2;

% Dt(1:168) = [0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18 ...
%  0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18 ...
%  0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18 ...
%  0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18 ...
%  0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18 ...
%  0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18 ...
%  0.18	0.18	0.18	0.18	0.18	0.66 ...
%  0.66	0.18	0.18	0.18	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18];
% 
% Dt = repmat(Dt',1,2);

%valCoeffBase = 0.5*N*max(Dt)/(min(Dt)*9.*sqrt(var(Dt)));%*2.5*mean(Dt)); % valor empirico asociado a la valoracion que dan los usuarios al recurso en una hora determinada.
%valCoeffBase = N*max(Dt)/(5*min(Dt));%*2.5*mean(Dt)); 
valCoeff = 20;%valCoeffBase;%-2; 110 w/inc, 60 w/o inc

%Minimum energy consumption for all users during a time lapse.

%q_min is a matriz of (N, T_+1)
q_min = repmat(((mean((min(Dt))/1000)./m)), 1,T_+1);%.12. Value in p.u, first converted to kW, then normalized with the pop mass


pt = Dt./max(Dt)*valCoeff; %TODO: Analize pt and effect of coefficients in population stable state.

% number of strategies
T_ = length(Dt); % should be any value and work normal as well.

% valuation parameters of all agents
% valoracion heterogenea entre poblaciones
alpha_ef = zeros(N,T_);
for i=1:N/2
    alpha_ef(i,:) = pt(1:T_,1);%*(1+.5*i/N*0) + 0.*rand(1.T_);
end

%pt = Dt./max(Dt)*(valCoeff+2); para cambiar la valoración por grupos de
%usuarios

for i=(N/2+1):N
    alpha_ef(i,:) = pt(1:T_,2);%*(1+.5*i/N*0) + 0.*rand(1.T_);
end

% parametros de la func. de costo agregado
beta_ef = 2;
b = 0;

% Time of the activation of either incentives or attacks IS NOT doing
% anything in the code
%time_on = 2;
%time_off = 4;


% number of pure strategies per population
n = T_+1;%(24*7)+1;%25 % MC. Energia . o potencia? que se distribuye entre los agentes. por cada poblacion.

% simulation parameters
time = 30;

% initial condition
pot=ones(N,T_+1)/(1*(T_+1));

%pot(:,1:T_) = repmat(Dt,N,1);%
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
%U = utility(x);

%figure(3);
%clf
%plot(1:1:T_, U(1,1:T_))
%tit = strcat('Utility per population, % addcost: ', num2str(addcost));
%title(tit)
%name = strcat('utility','%',num2str(addcost));
%savefig(3, name, 'compact');

%% sumas de potencia para graficar
sum_x = zeros(T_,1);
sum_x0 = zeros(T_,1);
for i=1:(T_)
    sum_x(i) = sum( x(:, i) );
    sum_x0(i) = sum( mp*x0(:,i));
end

figure(4); 
clf
plot(1:1:T_, sum_x(1:T_,1), 1:1:T_, sum_x0(1:T_,1))
titPower = strcat('Power allocation in the society [kW]');%, % addcost: ', num2str(addcost));
%title(titPower)
xlabel('Time [min]') % x-axis label
ylabel(titPower) % y-axis label
legend('Final allocation','Starting allocation')
nameP = strcat('power',netType,num2str(N),social_incentives);
%hold on
%    plot(1:1:T_, sum_x0(1:T_,1))
savefig(4, nameP, 'compact');


%% Resource allocation verification test
resources_alloc_pu = zeros(G.P,1);
resources_allocated = zeros(G.P,1);
for count = 1 : G.P
    resources_alloc_pu(count,1) = sum(x_n(count,1:T_));
    resources_allocated(count,1)= sum(x(count,1:T_));
end

resources_alloc_pu
total_resources_alloc = sum(resources_allocated);
disp(['Allocated total energy: ' num2str(total_resources_alloc) ' kWh']);
disp(['Available total energy: ' num2str(sum(m)) ' kWh']);
disp(['Energy use ratio: ' num2str(total_resources_alloc/sum(m))]);

figure(7)
iter= 1:N;
plot(iter,theta, iter, w, iter, s)
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
