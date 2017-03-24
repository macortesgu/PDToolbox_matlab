function [F, I] = fitness_user_finite(S,Q,Prefs,T_,IncActive,op)
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

power = Q/1000;
theta = Prefs;
%addcost es un porcentaje de costo agregado en forma de incentivos
%negativos (mayor costo,menor utilidad) para los agentes con un perfil de  
%consumo mayor.

%popul = size(Q,2);

%r_base = ones(1,(n-1))/(n-1); %debe sumar 1;
%r = r_base*mp;


sumQ = sum( sum(power) );
beta = unit_cost(sumQ);
%sumQi = sum(power);
ipow = 1+ repmat(devices.power,1,size(Q,2))/1000;

preferred_choice = eq((S-1), squeeze(theta(:,T_,:)));%    zeros(size(Q));
devices_off = eq(S,ones(size(S)));

%for p = 1 : popul
        
    %F(:,popul) = max(beta*(f*theta(:,T_,p)-1),0);
    %(1 + op(index)*(r(l) - q_t)/q_t)  
 %F(:,p) = (beta*f*preferred_choice(:,p).*(power(:,p)+1)) - beta*(power(:,p));
%F = (beta*f*preferred_choice.*(power + 1)) - beta*(power); 
alpha = value(beta,op);

if(IncActive == 1)
    I = incentives(ipow,beta,devices_off);
else
    I = zeros(size(ipow));
end

 F = (alpha.*preferred_choice) - beta*(power) + I;   
%end

