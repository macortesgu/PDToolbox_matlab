% test finite population

% run game with a finite population
% population games tool box
clear

global G

% number of populations
P = 1;

% number of agents
N = 200;

% number of pure strategies per population
n = 3;

% mass of the population
m = 1;

% initial condition
x0 = [0.2 .7 0.1]'; 

% simulation parameters
iterations = 10000;


% structure with the parameters of the game
G = struct('N', N, 'n', n, 'f', @fitness1, 'x0', x0, 'ode', 'ode113', 'time', iterations, 'eta', 0.02, 'revision_protocol', @proportional_imitation);

G.R = 1;

% verify data of the game
G = definition(G);

% have to include the combination of protocols!!
% have to include different tests for all the protocols.
G.revision_protocol = @proportional_imitation;
G.run_finite();
G.graph()
G.graph_evolution()
graph_utility
pause


G.revision_protocol = @comparison2average;
G.run_finite();
G.graph()
G.graph_evolution()
graph_utility
pause


G.revision_protocol = @pairwise_comparison;
G.run_finite();
G.graph()
G.graph_evolution()
graph_utility
pause


G.revision_protocol = @logit_choice;
G.run_finite();
G.graph()
G.graph_evolution()
graph_utility
pause