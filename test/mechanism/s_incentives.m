function [I, theta_step]= s_incentives(theta,alpha,epsilon,op2,s,s_step, energy_cost)
% pwc: Wether device is off or on.
% op: opinion for aech user respect to its acceptance of the social
% incentive. Initial opinion shuold be the susceptibility.

rng default;  % for reproducibility
sigma = binornd(1,op2);

incSignal = 0; %Signal generated with the incentive, to indicate if it is desired a reduction or an increase in consumption
%persistent theta_step_acum;

flexibility = 0.3;

theta_step = sigma.*((s_step*incSignal - s_step*(~incSignal))); %step change in preferences

%theta_step_acum = theta_step_acum_prev + theta_step;

delta_theta_inc = s - (theta + theta_step');

vartheta_inc = calc_vartheta(delta_theta_inc,epsilon);

s_reward_factor = 0.1; %factor considering a plus of shifting preferences, mking the new preference more desirable than the natural one. Accounting for users change'histeresis

I = alpha.*vartheta_inc - energy_cost*(1 - s_reward_factor);

%TODO falta hacer que el cambio en theta sea permanente y se acumule.
     
end