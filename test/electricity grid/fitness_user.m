function F = fitness_user(x, p)

global beta_ef alpha_ef mp

power = x;
T_ = 24;
F = zeros(T_+1,1);%F = zeros(T_+1,1);

index = p;
n = 24+1; %por qué 24+1?

for l = 1 : 24
    q_t = power(index, l);
    sum_q = sum( power(:, l) );
    alpha = alpha_ef(index, l);
    F( l ) = alpha / (1+q_t) - 2 * beta_ef * (sum_q) ;
end
F(n) = 0;


