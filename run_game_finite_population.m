function run_game_finite_population(name)
% RUN_GAME_FINITE_POPULATION Solves the difference equation of a populaiton game  
% 
% SYNOPSIS: RUN_GAME_FINITE_POPULATION(name)
% 
% INPUT name: Name of the structure that represents the game
% 
% The solution of the game, namely the evolution of strategies X in time,
% are attached to the game's structure 
% 
% REMARKS Its better to execute first <a href="matlab: help definition">definition</a> and run the game using G.run_finite(). 
%         This function uses the global variable 'G' to define the game
%
% SEE ALSO add_path, definition, run_game
%
% For more information see: <a href="https://github.com/carlobar/PDToolbox_matlab/">the GitHub's repository.</a>
% 
% Carlos Barreto, 04-11-16 
%
% Modified to run multipopulation games with finite number of agents, N.
% Mateo Cortés, March 2017.

global G

% load the structure of the game that calls the function
G = evalin('base', name);

% set initial strategy of each agent to satisfy the initial condition
protocol = func2str(G.revision_protocol);

if G.verb == true
    disp (['Running the ', protocol, ' revision protocol, period T=',num2str(G.period)]);
    tic
end

theta = squeeze(G.Prefs(:,G.period,:));
s = zeros(G.N, G.P);% randi(max(G.S),G.N,G.P); %%(max(G.S).*theta); %strategy matrix. Each col is a population.

% calculate the initial strategy of each agent given the proportions in x0


%for a = 1: G.P % for all populations. Current population is 'a'
%  h = 0;  
%    for i = 1: max(G.S(a))
%        p = floor(G.N * G.x0(a, i));
%        if ((p + h) <= G.N) && (p ~= 0)
%            s(h + 1: h + p, a) = i;
%            h = h + p;
%        end
        
%    end 
%end

% choose a random strategy to complete the initial state vector
%for a = 1: G.P
%    if h ~= G.N
%        s(h + 1: G.N, a) = unidrnd(G.S(a), 1, G.N - h);
%    end
%end

% set the number of iterations
t_max = floor(G.time);
    
T = 1:1:t_max;
%X = zeros( t_max, G.P*max(G.S));
QT = zeros(t_max, G.N*G.P);
Fitness = zeros(t_max, G.N*G.P);
Inc = zeros(t_max, G.N*G.P);

% Number of agents that update their strategy at each time instant, for
% each population

alarm = poissrnd( G.N * G.R, G.P, t_max);
%Q = zeros(G.N,G.P);

incentivablePeriods = Iper(); 


if(incentivablePeriods(G.period)==1)
    IncActive = 1;
else
    IncActive = 0;
end

%% Opinion dynamics

initial_state = linspace(0.8,0.2,G.P);% %zeros(N,1)initial beliefs of users
%c_eps = 1e-4;%1e-6; error aceptable para la convergencia
%opEps = 0.4; Bounded confidence parameter for HK op dyn model

% susceptibility (1-stubborness) theta of the agents
susceptibility = ((logspace(2, 3, G.P)/1000))';%column vector of susceptibility
susceptibility(G.P) = 0.35;
%define self confidence
w1 = logspace(0, 3, G.P)/1000;
w2 = fliplr(w1);
w = (w1+w2)/2;
w = w';

if(G.period == G.T_)
    figure(7)
    iter= 1:G.P;
    plot(iter,susceptibility, iter, w, iter, initial_state)
    legend('susceptibility','self-conf','i.c.','Location','NorthEastOutside')
    xlabel('Users');
    ylabel('Parameter value')
    % savefig(7, 'parameters', 'compact');
end

persistent Aw

if isempty(Aw)
   Aw = step_1_graph(G.P);
end

opinion = zeros(G.P,t_max);
opinion(:,1) = initial_state;

opinion2 = zeros(G.P,t_max); %opinion regarding susceptibility to change preferences due to social incentive.
opinion2(:,1) = susceptibility; %initial conditions for 2nd dynamic are the susceptibilites of the agents.

%opinion = bullo2(initial_state,P,Aw,theta,tMax,w, c_eps, 'plot');
%[opinion]=bullo2(s,N,A,theta,tMax,w, c_eps, varargin)

%theta_step = zeros(G.P,1);
s_step = 1/max(G.S);
%% Everything else  
for t = 1: t_max
    
    %Apagar los incentivos 
   % if(t>500)
   %     IncActive = 0;
   % end
    
    %% Opinion dynamics
    
    for i=1:G.P
       
        [~, J, Vneigh] = find(Aw(i,:)); % indexes and value of non-zero matriz elements.
        sumNeighbour = sum( opinion(J,t).*Vneigh');  %sum of neighbours opinion and confidence weights. Confidence is random
        sumNeighbour2 = sum( opinion2(J,t).*Vneigh');
        
        opinion(i, t+1)=(1-susceptibility(i))*opinion(i,1)+susceptibility(i)*w(i)*opinion(i,t)+susceptibility(i)*(1-w(i))*sumNeighbour;
        opinion2(i, t+1)=(1-susceptibility(i))*opinion2(i,1)+susceptibility(i)*w(i)*opinion2(i,t)+susceptibility(i)*(1-w(i))*sumNeighbour2;
        
    end
    % % Check if we have reached an equilibrium
    %opinionDiff = diff(X(:,t:t+1),1,2);
    %if (norm(opinionDiff, NORM_TYPE) < c_eps)
    %    disp(['[Bullo stubborness model] Reached equilibrium after ' num2str(t) ' rounds.']);
    %    break;
    %end
    
    
    
    %% Population Dynamics
    % update society state
    %x = state(s);
    Q = power_state(s,G.S'); %state(s);
    %for demand management problem, state of interest is
    %the power consumption associated to the users preferences of
    %wether to use or not their devices
    
    % find the current payoff of each strategy
    %F = zeros(G.S(1), a) ;
  %  if (t == 1)
 %       theta = squeeze(G.Prefs(:,G.period,:)); %2D preferences matrix for the current period
%addcost es un porcentaje de costo agregado en forma de incentivos
%negativos (mayor costo,menor utilidad) para los agentes con un perfil de  
%consumo mayor.
    %if t>1
    %    theta = max(theta + theta_step',0); %users' preferences are shifted due to the incentive, each iteration
    %end

 
    dualOpinion = [opinion(:,t) , opinion2(:,t)];
    
    s_norm = (s./(G.S'))-s_step;

    [F, I, theta_step, deltatheta]= G.f(s_norm, Q, theta, IncActive, dualOpinion,t,s_step); %f(x, a);
     if (t == 299)
         disp(deltatheta)
     end
    
    %update_agents = zeros(G.P,1);
    for a = 1: G.P
        % select users to update their actions at random
        if mod(t,5)==0
             update_agents = 1:G.N;
             limit = G.N;
        else
        update_agents = unidrnd(G.N, 1, alarm(a,t));
        limit = alarm(a,t);
        end
        
        % procedure to update the strategy of each agent
        s_update = s;
        for k=1 : limit
            i = update_agents(k);
       % t
        
            s_update(i,a) = G.revision_protocol(F, s, i, a, G.Prefs, IncActive,dualOpinion,deltatheta(i,a)); %(F(i,a), x, s, i, a)
        end
         s = s_update;   
    end
    

    
    %X(t, :) = x(:)';
    QT(t,:) = Q(:)';
    Fitness(t,:) = F(:)';
    Inc(t,:) = I(:)';
    
    
    
  
end

if G.verb == true
    toc
    disp(' ')
end
%% Opinion dynamcics
%opinion(:,t+2:end) = [];
%opinion(:,t_max+1) = []; 
%opinion2(:,t_max+1) = [];

    % Plot
    if(G.period == G.T_)
        plotOpinions2(opinion,21,'Users electricity valuation');
        hold on
        plotOpinions2(opinion2,22, 'Willingness to collaborate with DR program');
        hold off
    end

%% 
%G.X = X;
G.T = T;
G.Q = QT;
G.F = Fitness;
G.Inc = Inc;
G.opinion = opinion;
G.opinion2 = opinion2;

% save changes on the structure in the workspace
assignin('base', name, G);
clear G
