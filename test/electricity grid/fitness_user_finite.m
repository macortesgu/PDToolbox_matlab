function F = fitness_user_finite(S,Q,Prefs,T_)
% input:
% Q = population state (power consumption) (G.N,G.P)
% Prefs = 3D Matrix of preferences for users, of size (G.N, T_, G.P)
% T_ = number of periods for the secuence of games.
% output:
% F = payoff for each agent according to its strategy, for all agents in a
% period T_. Size (G.N,G.P)

%global beta_ef alpha_ef mp op

power = Q/1000;
theta = Prefs;
%n = 3; %debe coincidir con la definición de la variable en test_RA.m
F = zeros(size(Q)); % puede ser diferente para cada agente que juega una estrategia particular, para cada población.

%addcost es un porcentaje de costo agregado en forma de incentivos
%negativos (mayor costo,menor utilidad) para los agentes con un perfil de  
%consumo mayor.

popul = size(Q,2);

%r_base = ones(1,(n-1))/(n-1); %debe sumar 1;
%r = r_base*mp;
f = 2; %factor de preferencia, cuántas veces por encima del costo unitario valora cada persona una unidad de energia consumida

sumQ = sum( sum(power) );
beta = unit_cost(sumQ);


preferred_choice = eq((S-1), squeeze(theta(:,T_,:)));%    zeros(size(Q));


for p = 1 : popul
        
    %F(:,popul) = max(beta*(f*theta(:,T_,p)-1),0);
    %(1 + op(index)*(r(l) - q_t)/q_t)  
 F(:,p) = (beta*f*preferred_choice(:,p).*(power(:,p)+1)) - beta*(power(:,p));
%F = (beta*f*preferred_choice.*(power + 1)) - beta*(power);
    
end

