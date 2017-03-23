%% Step 2 opdyn

susceptibility = linspace(1,0.2,P);% %zeros(N,1)initial beliefs of users
tMax = 1e3; %Max. iterations for op dynamics
c_eps = 1e-4;%1e-6; error aceptable para la convergencia
%opEps = 0.4; Bounded confidence parameter for HK op dyn model

%% susceptibility (1-stubborness) theta of the agents
theta = (1-(logspace(0, 3, P)/1000))';%column vector of susceptibility
%define self confidence
w1 = logspace(0, 3, P)/1000;
w2 = fliplr(w1);
w = (w1+w2)/2;
w = w';
%fprintf(1,'å: %3.2f\n', opEps);   
%[~, opinion] = hkLocal(A, s, opEps, tMax, c_eps, 'plot');

opinion = bullo2(susceptibility,P,Aw,theta,tMax,w, c_eps, 'plot');

social_incentives = 'sInc';
switch social_incentives
    case 'sInc'
        op = opinion(:,end); %zeros(1,T_);%
        
    otherwise
        op = zeros(P,1); %No use of social incentives
        social_incentives = 'noneInc';
end



%save('Caso1.mat','Gw','opinion'); 