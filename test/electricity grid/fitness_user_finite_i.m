function F = fitness_user_finite_i(preferred_choice,beta,Q,i,p)

% input:
% Q = population state (power consumption) (G.N,G.P)
% Prefs = 3D Matrix of preferences for users, of size (G.N, T_, G.P)
% T_ = number of periods for the secuence of games.
% output:
% F = payoff for each agent according to its strategy, for all agents in a
% period T_. Size (G.N,G.P)

global devices
%r_base = ones(1,(n-1))/(n-1); %debe sumar 1;
%r = r_base*mp;
ipow = 1 + devices.power(i)/1000;
alpha = value(beta);%beta*f;
alpha = alpha(p);
device_off = eq(Q,0);
%F = max(beta*(f*theta - 1),0);
I = incentives(ipow,beta,device_off);
%F = (alpha*preferred_choice*((Q/1000)+1)) - beta*(Q/1000);
F = (alpha*preferred_choice) - beta*(Q/1000) + I;
    %(1 + op(index)*(r(l) - q_t)/q_t)  

