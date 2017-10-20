function s_i = logit_choice(F, s, i, p, Prefs, IncActive,dualOp,deltatheta, s_step)%(F, z, s, i, p)
% LOGIT_CHOICE Computes the differece equation that describes the update of
%         the populations' state following the logit choice revision protocol. 
%         This revision protocol leads to the logit dynamics with a large 
%         number of agents
% 
% SYNOPSIS: s_i = LOGIT_CHOICE(F, z, s, i, p)
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
%          pairwise_comparison, proportional_imitation
%
% For more information see: <a href="https://github.com/carlobar/PDToolbox_matlab/">the GitHub's repository.</a>
% 
% Carlos Barreto, 04-11-16 
% MC mod for pairwise logit protocol.

global G

% j = unidrnd( G.S(p) );
% 
% 
% eta = G.eta;
% F_ = exp( F(1:G.S(p) ) / eta );
% F_ref = sum(F_(:));
% 
% rho_ij =  F_(j) / F_ref;
% 
% % prob generator
% change = ceil(rand - 1 + rho_ij / G.R);
% 
% if change == 1
% 	s_i = j;
% else
% 	s_i = s(i);
% end
% 
% if j == 3
% 	h=1;
% end
%% Replacing part


%s_inorm = s(i,p)/max(G.S);
theta = Prefs(i,G.period,p);

num_s = max(G.S); % as 1st strategy means zero consumption;

%j = unidrnd( num_s );

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


s_try = s;
s_try(i,p) = j;
Q_try = power_state(s_try,G.S');
sumQtry = sum(sum(Q_try));

s_try_norm = s_try(i,p)/G.S(p) - 1/max(G.S);

i_norm = s(i,p)/G.S(p) - 1/max(G.S);


pi_j = fitness_user_finite_i(theta,s_try_norm,i,sumQtry,IncActive,dualOp(p,:), s_step);%fitness_user_finite_i(preferred_choice,beta,Qj); op(p) que toma si recibe dualOpinion?

pi_i = fitness_user_finite_i(theta,i_norm,i,sumQtry,IncActive,dualOp(p,:), s_step);

eta = G.eta;
F_ = exp( (pi_j/1000) / eta );
F_ref = F_ + exp((pi_i/1000) / eta);

rho_ij =  F_ / F_ref;


% prob generator
change = ceil(rand - 1 + rho_ij / G.R);

% if(rho_ij==0)
%     change = 0;
% else
%     change = binornd(1,rho_ij);%%1 & rho_ij;ceil(rand - 1 + rho_ij / G.R);
% end

%if(deltatheta==0 && change==1)
%    warning('Agente en estado óptimo pero cambia estrategia')
%end

if change == 1 %==1
	s_i = j;
else
	s_i = s(i,p);
end
