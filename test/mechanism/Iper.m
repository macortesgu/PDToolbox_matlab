function Iperiods = Iper(Iper)

persistent IncPeriods
% Parse input
switch nargin
    case 1
       IncPeriods = Iper; 
       Iperiods = 'done';
    otherwise
       Iperiods = IncPeriods;
end
  
end 