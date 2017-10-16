function I = incentives(ipow,beta,devices_off)
% I = incentives(ipow,beta,devices_off,validPeriod)
% ipow: Power reduction due to user specific device' preference change. It
% is equal to the corresponding device nominal power usage.
% beta: Unit cost according to the total energy consumption
% devices_off: Boolean to restrict the incentive to the case the user's devices
% are turned off
% validPeriod: Boolean vector indicating if the incentives are given for
% each period (divisions for a week).

gamma = 5*beta; %factor of additional surplus for changinf preference

%validPeriod = Iper();%period; %Analyze if consumption for the period is below or above average. If it is below average, then
% no incentive is needed as it would not contribute to flatten the weekly
% profile.

I = gamma.*(ipow).*devices_off;