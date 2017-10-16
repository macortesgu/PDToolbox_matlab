function vartheta = calc_vartheta(deltatheta,epsilon)
% vartheta = calc_vartheta(deltatheta,epsilon)
% Calculates degree of fullfillment of user actual preferences, according
% to their strategies being played.
vartheta = exp(-abs(deltatheta)/epsilon);
%max(1 - abs(deltatheta)/epsilon,0);