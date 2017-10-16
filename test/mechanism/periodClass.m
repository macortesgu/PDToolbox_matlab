function label = periodClass(T_)

base = 'Periods';

switch T_
    case 168
        label = strcat(base,' [h]');
    case 10080
        label = strcat(base,' [min]');
    otherwise
        units = 10080/T_;
        label = strcat(base, ' [',num2str(units),' min]');
end