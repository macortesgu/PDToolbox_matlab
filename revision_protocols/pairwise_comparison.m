function s_i = pairwise_comparison(F, s, i, p, Prefs, IncActive,op,deltatheta) % (F, z, s, i, p)
% PAIRWISE_COMPARISON Computes the differece equation that describes the update
%         of the populations' state following the pairwise comparison revision
%         protocol. This revision protocol leads to the Smith dynamics with
%         a large number of agents
% 
% SYNOPSIS: S_I = PAIRWISE_COMPARISON(F, z, s, i, p)
% 
% INPUT F: Vector with the payoff of each strategy
%       z: Vector with the society's state (distribution of strategies)
%       s: Vector with the strategy of each agent
%       i: ID of the agent that makes an update of its strategy
%       p: Population ID. The current version supports finite games with
%          more than one population
% 
% OUTPUT s_i: Update of the agent's strategy
%
% SEE ALSO definition, run_game_finite_population, comparison2average, 
%          logit_choice, proportional_imitation
%
% For more information see: <a href="https://github.com/carlobar/PDToolbox_matlab/">the GitHub's repository.</a>
% 
% Carlos Barreto, 04-11-16 

global G 

%s_inorm = s(i,p)/max(G.S);
theta = Prefs(i,G.period,p);

num_s = max(G.S)+1; % as 1st strategy means zero consumption;



if deltatheta > 0
   % deltatheta
   % [1 s(i,p)]
    j = randi([1 round(s(i,p))],1);%unidrnd( G.S(p) );    
else
    if (deltatheta==0)
        j = randi(num_s,1); %if deltatheta is 0, means agent has best response strategy, so there are no limitations on random strategy selection
        
    else
      %  deltatheta
       % [s(i,p) num_s]
    j = randi([round(s(i,p)) round(num_s)],1);
    end
end

%pi_i = F(i,p);

s_try = s;
s_try(i,p) = j;
Q_try = power_state(s_try,G.S');
sumQtry = sum(sum(Q_try));

s_try_norm = s_try(i,p)/G.S(p) - 1/max(G.S);

i_norm = s(i,p)/G.S(p) - 1/max(G.S);


pi_j = fitness_user_finite_i(theta,s_try_norm,i,sumQtry,IncActive,op(p));%fitness_user_finite_i(preferred_choice,beta,Qj);

pi_i = fitness_user_finite_i(theta,i_norm,i,sumQtry,IncActive,op(p));

rho_ij = min(max(pi_j - pi_i, 0)/200,1);%max(pi_j - pi_i, 0);

% prob generator

if(rho_ij==0)
    change = 0;
else
    change = binornd(1,rho_ij);%%1 & rho_ij;ceil(rand - 1 + rho_ij / G.R);
end

if(deltatheta==0 && change==1)
    warning('Agente en estado óptimo pero cambia estrategia')
    disp(pi_j)
    disp(pi_i)
    disp(j)
    disp(s(i,p))
    disp(s_try_norm)
    disp(theta)
end

if change == 1 %==1
	s_i = j;
else
	s_i = s(i,p);
end

