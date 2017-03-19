function F = fitness_user_finite_i(preferred_choice,beta,Q)

% input:
% Q = population state (power consumption) (G.N,G.P)
% Prefs = 3D Matrix of preferences for users, of size (G.N, T_, G.P)
% T_ = number of periods for the secuence of games.
% output:
% F = payoff for each agent according to its strategy, for all agents in a
% period T_. Size (G.N,G.P)


%r_base = ones(1,(n-1))/(n-1); %debe sumar 1;
%r = r_base*mp;
f = 2; %factor de preferencia, cuántas veces por encima del costo unitario valora cada persona una unidad de energia consumida

%F = max(beta*(f*theta - 1),0);
F = (beta*f*preferred_choice*((Q/1000)+1)) - beta*(Q/1000);
    %(1 + op(index)*(r(l) - q_t)/q_t)  

