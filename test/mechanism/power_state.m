function Q = power_state(s,Stot)
% POWER_STATE Computes the power consumption for each user(population of devices) during the period.
% 
% SYNOPSIS: z = STATE(s)
% 
% INPUT s: Vector with the strategy of each agent. //Matrix with the
% strategy of each agent for all populations. SIZE (G.N x G.P)
%       
% OUTPUT z: Vector with the social state //Matrix with the social state
%
% SEE ALSO definition, run_game_finite_population, comparison2average, 
%          logit_choice, pairwise_comparison, proportional_imitation
%


global devices

% find the current social state 
%z = zeros(max(G.S), G.P);

%devicesState = (s-1) & ones(size(s));
pow = repmat(devices.power,1,size(s,2));
Q =pow.*(s./Stot);


