function alpha = value(beta)
global op
f = 5; %factor de preferencia, cuántas veces por encima del costo unitario valora cada persona una unidad de energia consumida
alpha = f*beta*op';% ones(1,length(op))