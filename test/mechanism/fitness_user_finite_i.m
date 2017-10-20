function F = fitness_user_finite_i(theta, s_try_norm, i,sumQtry,IncActive,op, s_step)

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

epsilon = 0.15;
minq = 0.01;

deltatheta = s_try_norm - theta;
vartheta = calc_vartheta(deltatheta,epsilon);
device_off = s_try_norm<minq;

beta = unit_cost(sumQtry);
ipow = 1 + devices.power(i)/1000;
alpha = value(beta,op(1,1));%beta*f;

energy_cost = beta*(Qj/1000);

if(IncActive == 1)
    [Isoc, theta_step] = s_incentives(theta,alpha,epsilon,op(1,2),s_try_norm,s_step,energy_cost);
    I =Isoc; %+ incentives(ipow,beta,device_off) ;
else
    %theta_step = 0;
    I = 0;%zeros(size(ipow));
end

F = (alpha*vartheta) - energy_cost + I;
    

