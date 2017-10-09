function F = fitness_user_finite_i(theta, s_try_norm, i,sumQtry,IncActive,op)

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


Qj = devices.power(i)*(s_try_norm);%Power consumption for device i wether it is on or off

epsilon = 0.2;
minq = 0.04;

deltatheta = s_try_norm - theta;
vartheta = calc_vartheta(deltatheta,epsilon);
device_off = s_try_norm<minq;

beta = unit_cost(sumQtry);
ipow = 1 + devices.power(i)/1000;
alpha = value(beta,op);%beta*f;

if(IncActive == 1)
    I = incentives(ipow,beta,device_off);
else
    I = 0;%zeros(size(ipow));
end

F = (alpha*vartheta) - beta*(Qj/1000) + I;
    

