function F = fitness_user(x, p)

global beta_ef alpha_ef mp addcost T_ q_min op

power = x;
F = zeros(T_+1,1);

%addcost es un porcentaje de costo agregado en forma de incentivos
%negativos (mayor costo,menor utilidad) para los agentes con un perfil de  
%consumo mayor.

addcostPU = addcost/100; % normalizar

index = p;
n = T_+1; %debe coincidir con la definición de la variable en test_RA.m

r_base = ones(1,T_)/T_; %debe sumar 1;
r = r_base*mp;

for l = 1 : T_
    q_t = power(index, l);
    sum_q = sum( power(:, l) );
    alpha = alpha_ef(index, l);
    %Constraints check
   % if q_t >= (q_min+0.03)
        F(l) = (alpha/(1+q_t))*(1 + op(index)*(r(l) - q_t)/q_t) - 2*beta_ef*(sum_q);
      %  F(l) = (alpha/(1+q_t))*(1 + op(index)*(r(l) - q_t)/q_t) - (1+addcostPU*alpha)*beta_ef*(sum_q+q_t);
  %  else
  %      F(l) = 0;
  %  end
end


F(n) = 0;


