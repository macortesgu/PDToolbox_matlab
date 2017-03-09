function F = fitness_user(x, p)

global beta_ef alpha_ef mp addcost T_

power = x;
%T_ = 24;
F = zeros(T_+1,1);

%addcost es un porcentaje de costo agregado en forma de incentivos
%negativos (mayor costo,menor utilidad) para los agentes con un perfil de  
%consumo mayor.

addcostPU = addcost/100; % normalizar

index = p;
n = T_+1; %debe coincidir con la definición de la variable en test_RA.m

for l = 1 : T_
    q_t = power(index, l);
    sum_q = sum( power(:, l) );
    alpha = alpha_ef(index, l);
    F(l) = alpha/(1+q_t) - (1+addcostPU*alpha)*2*beta_ef*(sum_q) ;
    %F( l ) = alpha / (1+q_t) - 2 * beta_ef * (sum_q) ;
end
F(n) = 0;


