function p = unit_cost(Q)
% Unitary price function.
%INPUT = total energy Q in [W], consumed by a society.
%OUTPUT = unitary price per kWh consumed.
% the function assumes equal price for different time periods, and a
% dependency only on the quantity of energy generated.

pH = 300;
pC = 350;
pG = 500;
Hmax = 30;
Cmax = 50;
Gmax = 60;

QkWh = Q/1000;

if QkWh <= Hmax
    p = pH;
elseif (QkWh>Hmax && QkWh<=Cmax)
    p = pC;
elseif (QkWh>Cmax && QkWh<=Gmax)
    p = pG;
else
    p = pG*1.1;
end