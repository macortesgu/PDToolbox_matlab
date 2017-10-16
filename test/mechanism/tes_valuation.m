path(path, '../../revision_protocols')
path(path, '../../graphs')
path(path, '../../dynamics')
path(path, '../../')


clear

%l = 10;
ith = 1;


global G

global beta_ef alpha_ef mp N T_ addcost

N = 10; %Number of users

size_Pt = 1:2:150; %sweep span
addcost = 0;
%ptCoeff = 4; % valor empirico asociado a la valoracion que dan los usuarios al recurso en una hora determinada.

% Definition of the electricity variables
Dt = [0.18	0.18	0.18	0.18	0.18	0.66 ...
 0.66	0.18	0.18	0.18	0.30	0.30 ...
 0.20	0.31	0.36	0.18	0.28	0.28 ...
 0.47	0.45	0.39	0.28	0.18	0.18];

% Dt = [0.18	0.18	0.18	0.18	0.18	0.8 ...
%  0.66	0.18	0.18	0.4	0.30	0.30 ...
%  0.20	0.31	0.36	0.18	0.28	0.28 ...
%  0.47	0.45	0.39	0.28	0.18	0.18];

T_ = length(Dt);

%Resulting values of interest after the parametric sweep
P_max = zeros(T_, length(size_Pt));
P_min = zeros(T_, length(size_Pt));
P_l = zeros(T_, length(size_Pt));

U_l = zeros(T_, length(size_Pt));


    
    for k=1 : length(size_Pt)
        
        disp (['Running  Valuation=', num2str(size_Pt(k)), '... ']);
        
        ptCoeff = size_Pt(k); %sweep
        
        pt = Dt./max(Dt)*ptCoeff; %TODO: Analize pt and effect of coefficients in population stable state.
        
        % valuation parameters of all agents
        alpha_ef = zeros(N,T_);
        for i=1:N
            alpha_ef(i,:) = pt(1:T_);%*(1+.5*i/N*0) + 0.*rand(1.T_);
        end
        
        % parametros de la func. de costo agregado
        beta_ef = 1;
        
        
        
        % definition of the game structure
        
        % number of populations
        P = N;
        
        % number of pure strategies per population
        n = 25;
        mp = 6.8825; %14
        m = ones(P, 1) * mp;
        dyn = {'rd'};
        
        % simulation parameters
        time = 30;
        
        %Initial condition
        %pot_r = ones(N,T_+1)*mp/(T_+1);
        %x0 = pot_r/mp;
        pot = ones(N,T_+1)/(T_+1);
        x0 = pot;
        
        % structure with the parameters of the game
        G = struct('P', P, 'n', n, 'f', @fitness_user_inefficient, 'ode', 'ode45', 'time', time, 'tol', 0.00001, 'x0', x0, 'm', m);
        
        % random initial condition
        %G = struct('P', P, 'n', n, 'f', @fitness_user, 'ode', 'ode23s', 'time', time, 'step', 0.00001);
        
        % verify data of the game
        G = definition(G);
        
        G.step = .01;
        
        
        % run game
        G.dynamics = dyn;
        G.run()
        T_dyn = G.T;
        X_dyn = G.X;
        
        
        
        % extract matrix of strategies
        %n = max(G.S);
        x_n = vec2mat(X_dyn(end, :), n);
        x = zeros(G.P, n);
        
        for p = 1 : G.P
            x(p, :) = x_n(p, :) * G.m(p);
        end
        
        %U_i = utility_incentives(x_i);
        U = utility(x);
        
        %Save step information in results variables
        for l=1:T_
        
        U_l(l,k) = sum( U(:, l) );
        
        P_l(l,k) = sum( x(:, l) );
        
        end
        
        %P_max = zeros(1, length(size_Pt));
        %P_min = zeros(1, length(size_Pt));
        
        
        
        
        
        
        % run the simulations with the inneficient case
        %G.f = @fitness_user_inefficient;
        %G.run()
        %X_dyn = G.X;
        
        % extract matrix of strategies
        %n = max(G.S);
        %x_n = vec2mat(X_dyn(end, :), n);
        %x = zeros(G.P, n);
        
        %for p = 1 : G.P
        %    x(p, :) = x_n(p, :) * G.m(p);
        %end
        
        %U = utility(x);
        
        %U_nash(k) = sum( sum( U(:, l) ) );
        %X_nash(k) = sum( sum( x(:, l) ) );
        
        %figure(3); plot(1:1:24, U(ith, 1:24), 1:1:24, U_i(ith, 1:24), 'r')
        %figure(4); plot(1:1:24, x(ith, 1:24), 1:1:24, x_i(ith, 1:24), 'r')
        
    end

figure(1)
clf
hold on
for l=1:T_
    plot(size_Pt, U_l(l,:) )
end
title('Utilidad en función de la valoración del recurso');
hold off


figure(2)
clf
hold on
for l=1:T_
    plot(size_Pt, P_l(l,:) )
end
title('Potencia consumida en función de la valoración del recurso');
%plot(size_Pt, (size_Pt+1)./(2*size_Pt), '--k')
hold off













