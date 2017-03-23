function I = incentives(ipow,beta,devices_off)
% I = incentives(f,ipow)
gamma = 0.8*beta; %factor of additional surplus for changinf preference

I = gamma.*(ipow).*devices_off;