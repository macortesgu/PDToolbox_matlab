function alpha = value(beta,op)
%global op
f = 7; %factor de preferencia, cuántas veces por encima del costo unitario valora cada persona una unidad de energia consumida
alpha = f*beta*op';% ones(1,length(op))