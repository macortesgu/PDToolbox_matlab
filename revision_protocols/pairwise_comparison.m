function s_i = pairwise_comparison(F, s, i, p, Prefs) % (F, z, s, i, p)
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
%       p: Population ID. The current version only supports finite games with
%          one population
% 
% OUTPUT s_i: Update of the agent's strategy
%
% SEE ALSO definition, run_game_finite_population, comparison2average, 
%          logit_choice, proportional_imitation
%
% For more information see: <a href="https://github.com/carlobar/PDToolbox_matlab/">the GitHub's repository.</a>
% 
% Carlos Barreto, 04-11-16 

global G devices

%j = unidrnd( G.S(p) );

if s(i,p) == 1 %aparato apagado
    j = 2;
else 
    if (s(i,p) == 2)%aparato prendido
        j = 1;
    end
end

pi_i = F(i,p);

s_try = s;
s_try(i,p) = j;
Q_try = power_state(s_try);
sumQtry = sum(sum(Q_try));
beta = unit_cost(sumQtry);
theta = Prefs(i,G.period,p);
Qj = devices.power(i)*(s_try(i,p)-1);%Power consumption for device i wether it is on or off

preferred_choice = eq(theta,s_try(i,p)-1);

pi_j = fitness_user_finite_i(preferred_choice,beta,Qj,i,p);%fitness_user_finite_i(preferred_choice,beta,Qj);

rho_ij = max(pi_j - pi_i, 0)/1000;

% prob generator
change = 1 & rho_ij;%ceil(rand - 1 + rho_ij / G.R);

if change >= 1 %==1
	s_i = j;
else
	s_i = s(i,p);
end

