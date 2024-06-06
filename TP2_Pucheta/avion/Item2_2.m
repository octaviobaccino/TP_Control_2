clear; close all; clc

%--------------------------------------------------------------------------
% Alumno: Octavio Baccino
% Profesor: Julian Pucheta
% Trabajo Práctico 2
%--------------------------------------------------------------------------
% Caso de estudio 2. Sistema de cuatro variables de estado
%--------------------------------------------------------------------------
% Ítem [2]
% -Proponer un controlador en tiempo discreto en variables de estado para que 
% el proceso evolucione en los rangos de validez del modelo, es decir donde 
% los ángulos y el valor de la acción de control en valor absoluto son menores 
% a la unidad. -Asumiendo que no puede medirse el ángulo ?, pero sí el 
% ángulo ? y la altura, proponer un esquema que permita lograr el objetivo de 
% control. 
% -Establecer el valor del tiempo de muestreo más adecuado para 
% implementar el diseño en un sistema micro controlado. -Determinar el efecto 
% de la nolinealidad en la acción de control, descripta en la Fig. 4, y 
% verificar cuál es el máximo valor admisible de la nolinealidad.  
%--------------------------------------------------------------------------

%% Controlador en tiempo discreto

% x = [alpha; phi; phi_p; h]

config.a = 0.07;
config.b = 5;
config.w = 9;
config.c = 150;

a = config.a;
b = config.b;
w = config.w;
c = config.c;

Ac = [-a a 0 0; % Matriz de estados
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];

Bc = [0; 0; w^2*b; 0]; % Matriz de entradas 

Cc = [0 0 0 1]; % Matriz de salidas

Dc = 0; % Matriz de transmision directa

Q = diag([1e1 1e0 1e1 1e0 1e-2]);

R = 1e3;

Aac =  [Ac zeros(4,1); -Cc 0];

Bac = [Bc; 0];

% Discretizacion del sistema

Ts = 0.01;

sys_c = ss(Ac,Bc,Cc,Dc);
sys_d = c2d(sys_c,Ts,'zoh');

A = sys_d.A; % A discreta

B = sys_d.B; % B discreta

C = sys_d.C; % C discreta

D = sys_d.D; % D discreta

% Matrices ampliadas

Aa = [A zeros(4,1); -C*A 1]; % A discreta ampliada

Ba = [B; -C*B]; % B discreta amplaida

Q = diag([1e0 1e1 1e3 1e-2 1e-8]);

R = 1e-5;

[Ka,S,P] = dlqr(Aa,Ba,Q,R); % LQR en tiempo discreto

% K = matriz de ganacias
% S = solucion a la ecuacion de Riccati asociada
% P = polos de lazo cerrado del sistema

% G=inv(C*inv(eye(4)-A+B*K)*B);

fprintf('Polos de lazo cerrado:');
P
Ka

%% Simulacion

ref = 100;

x = [0; 0; 0; -500];

tsim = 100;

t_v = 0:Ts:tsim-Ts;

alpha = zeros(1,round(tsim/Ts));
phi = zeros(1,round(tsim/Ts));
phi_p = zeros(1,round(tsim/Ts));
h = zeros(1,round(tsim/Ts));

u = zeros(1,round(tsim/Ts));
ve = 0;
ve_k = 0;

for i=1:round(tsim/Ts)
    
    Y = C*x;
    y(i) = Y;
    ve = ref - y(i); %Ref no nula sólo para el desplazamiento
    ve_k = ve_k+ve;
    
    u(i)= -Ka(1:4)*x - Ka(5)*ve_k;

    x = modelo_avion(Ts,x,u(i),config);
    
    alpha(i)=x(1);
    phi(i)=x(2);
    phi_p(i)=x(3);
    h(i)=x(4);
    
end


%% Plots

figure(1);

subplot(3,2,1);
plot(t_v,alpha);
grid on; 
title('\alpha');
hold on;

subplot(3,2,2);
plot(t_v,phi);
grid on;
title('\phi');
hold on;

subplot(3,2,3); 
plot(t_v,phi_p);
grid on;
title('\phi_p');
hold on;

subplot(3,2,4);
plot(t_v,h);
grid on;
title('Altura h');
hold on;

subplot(3,1,3);
plot(t_v,u);
grid on;
title('Acción de control');
xlabel('Tiempo [s]');


