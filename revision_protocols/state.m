function z = state(s)
% STATE Computes the strategy's proportion (social state)
% 
% SYNOPSIS: z = STATE(s)
% 
% INPUT s: Vector with the strategy of each agent. //Matrix with the
% strategy of each agent for all populations.
%       
% OUTPUT z: Vector with the social state //Matrix with the social state
%
% SEE ALSO definition, run_game_finite_population, comparison2average, 
%          logit_choice, pairwise_comparison, proportional_imitation
%
% For more information see: <a href="GitHub: web('https://github.com/carlobar/PDToolbox_matlab/')">the GitHub's repository.</a>
% 
% Carlos Barreto, 04-11-16 


global G

% find the current social state 
z = zeros(max(G.S), G.P);

for a = 1:G.P
    for i = 1: max(G.S)
        z(i,a) = sum(s(:,a) == i) / G.N;
    end
end
