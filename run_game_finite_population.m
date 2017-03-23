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

s = ones(G.N, G.P); %strategy matrix. Each col is a population. Ones means all agents off.

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

% Number of agents that update their strategy at each time instant, for
% each population

alarm = poissrnd( G.N * G.R, G.P, t_max);
%Q = zeros(G.N,G.P);


   
for t = 1: t_max
    
    % update society state
    %x = state(s);
    Q = power_state(s); %state(s);
    %for demand management problem, state of interest is
    %the power consumption associated to the users preferences of
    %wether to use or not their devices
    
    % find the current payoff of each strategy
    %F = zeros(G.S(1), a) ;
    F = G.f(s, Q, G.Prefs, G.period); %f(x, a);
    
    %update_agents = zeros(G.P,1);
    for a = 1: G.P
        % select users to update their actions at random
        update_agents = unidrnd(G.N, 1, alarm(a,t));
        
        % procedure to update the strategy of each agent
        s_update = s;
        for k=1 : alarm(a,t)
            i = update_agents(k);

            s_update(i,a) = G.revision_protocol(F, s, i, a, G.Prefs); %(F(i,a), x, s, i, a)
        end
         s = s_update;   
    end
    

    
    %X(t, :) = x(:)';
    QT(t,:) = Q(:)';
    Fitness(t,:) = F(:)';
    
    
    
    
  
end

if G.verb == true
    toc
    disp([' '])
end


%G.X = X;
G.T = T;
G.Q = QT;
G.F = Fitness;

% save changes on the structure in the workspace
assignin('base', name, G);
clear G
