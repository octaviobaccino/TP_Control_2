function [X] = modelo_avion(Ts,estado,accion,config)

a = config.a;
b = config.b;
w = config.w;
c = config.c;

T0 = 2e-4; % debe ser menor que Ts

A = [-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];

B = [0; 0; w^2*b; 0];

C = [0 0 0 1];

D = 0;

X = [estado(1) estado(2) estado(3) estado(4)]';

alpha = zeros(1,round(Ts/T0));
phi = zeros(1,round(Ts/T0));
phi_p = zeros(1,round(Ts/T0));
h = zeros(1,round(Ts/T0));

for i=1:round(Ts/T0)
    
    alpha(i) = X(1);
    phi(i) = X(2);
    phi_p(i) = X(3);
    h(i) = X(4);
    
    Xp = A*X + B*accion;
    X = X + Xp.*T0;
    
end

end

