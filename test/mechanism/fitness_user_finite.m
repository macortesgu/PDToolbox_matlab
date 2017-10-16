function [F, I, theta_step, deltatheta] = fitness_user_finite(Snorm,Q,theta,IncActive,op,t,s_step)
% F (size (G.N,G.P) = fitness_user_finite(S,Q,Prefs,T_,IncActive)
% input:
% Q = population state (power consumption) (G.N,G.P)
% Prefs = 3D Matrix of preferences for users, of size (G.N, T_, G.P)
% T_ = number of periods for the secuence of games.
% output:
% F = payoff for each agent according to its strategy, for all agents in a
% period T_. Size (G.N,G.P)

%global beta_ef alpha_ef mp op
global devices

minq = 0.04;
power = Q/1000; %kWh
epsilon = 0.2; %rango de variación en la cual el usuario obtiene un beneficio, aunque no sea su preferencia exacta.

%popul = size(Q,2);

%r_base = ones(1,(n-1))/(n-1); %debe sumar 1;
%r = r_base*mp;


sumQ = sum( sum(power) );
beta = unit_cost(sumQ);
%sumQi = sum(power);
ipow = 1+ repmat(devices.power,1,size(Q,2))/1000;

%preferred_choice has into account the deviation from the desired load
%profile, being proportional to the difference between actual and desired
%preferences, if the difference is lesser than an elasticity limit. If
%outside the limit, it is zero.
%preferred_choice = eq((S-1), squeeze(theta(:,T_,:)));%    zeros(size(Q));


deltatheta = Snorm-theta;
vartheta = calc_vartheta(deltatheta,epsilon);

%if(t == 149)
 %   disp(deltatheta);
%end

devices_off = Snorm<minq;%eq(S,ones(size(S))); minq is the minimum proportion of perdio consumption considered as zero, in order to grant incentives

alpha = value(beta,op(:,1));

%pwc = rand(size(Q));
energy_cost = beta*(power);

if(IncActive == 1)
    [Isoc, theta_step] = s_incentives(theta,alpha,epsilon,op(:,2),Snorm,s_step,energy_cost);
    I = Isoc;%incentives(ipow,beta,devices_off) + ;
else
    theta_step = 0;
    I = zeros(size(ipow));
end

 F = (alpha.*vartheta) - energy_cost;% + I;   
%end


