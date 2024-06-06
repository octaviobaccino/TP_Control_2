clear; close all; clc

%--------------------------------------------------------------------------
% Alumno: Octavio Baccino
% Profesor: Julian Pucheta
% Trabajo Práctico 2
%--------------------------------------------------------------------------
% Caso de estudio 2. Sistema de cuatro variables de estado
%--------------------------------------------------------------------------
% Ítem [2]
% -Asumiendo que no puede medirse el ángulo ?, pero sí el ángulo ? y la 
% altura, proponer un esquema que permita lograr el objetivo de control.
% -Establecer el valor del tiempo de muestreo más adecuado para implementar 
% el diseño en un sistema micro controlado.
% -Determinar el efecto de la nolinealidad en la acción de control, descripta 
% en la Fig. 4, y verificar cuál es el máximo valor admisible de la nolinealidad. 
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

Cc = [0 1 0 0; 0 0 0 1]; % Matriz de salidas

Dc = 0; % Matriz de transmision directa

% Discretizacion del sistema

Ts = 0.1;
% Ts = 0.05;

sys_c = ss(Ac,Bc,Cc,Dc);
sys_d = c2d(sys_c,Ts,'zoh');

A = sys_d.A; % A discreta

B = sys_d.B; % B discreta

C = sys_d.C; % C discreta

D = sys_d.D; % D discreta

% Matrices ampliadas

Aa = [A zeros(4,1); -C(2,:)*A 1]; % A discreta ampliada

Ba = [B; -C(2,:)*B]; % B discreta amplaida

Q = diag([1e-3 1e8 1e8 1e4 1e-2*30]);

R = 1e-0;

[Ka,S,P] = dlqr(Aa,Ba,Q,R); % LQR en tiempo discreto

% K = matriz de ganacias
% S = solucion a la ecuacion de Riccati asociada
% P = polos de lazo cerrado del sistema

% G=inv(C(2,:)*inv(eye(4)-A+B*K)*B);

fprintf('Polos de lazo cerrado:');
P
fprintf('Matriz K ampliada:');
Ka

%% Observador

Ao = A'; % A del sistema dual
Bo = C'; % B del sistema dual
Co = B'; % C del sistema dual
Do = D'; % D del sistema dual

Qo = diag([1e-7 1e-7 1e-7 1e0]);

Ro = diag([1e-6 1e-6]);

[Ko,So,Po] = dlqr(Ao,Bo,Qo,Ro);

fprintf('Polos de lazo cerrado del observador:');
Po

%% Simulacion

x = [0; 0; 0; -500];
xobs = [0; 0; 0; 0];

tsim = 100;

t_v = 0:Ts:tsim-Ts;

ref = 100;

alpha = zeros(1,round(tsim/Ts));
phi = zeros(1,round(tsim/Ts));
phi_p = zeros(1,round(tsim/Ts));
h = zeros(1,round(tsim/Ts));

alpha_obs = zeros(1,round(tsim/Ts));
phi_obs = zeros(1,round(tsim/Ts));
phi_p_obs = zeros(1,round(tsim/Ts));
h_obs = zeros(1,round(tsim/Ts));

u = zeros(1,round(tsim/Ts));
acc = zeros(1,round(tsim/Ts));

ve = 0;
ve_k = 0;

Yobs = [0 ;0]; % salida inicial

for i=1:round(tsim/Ts)
    
    Y = C*x; % salida real
    y(i) = Y(2,1);
    ve = ref - y(i); %Ref no nula sólo para el desplazamiento
    ve_k = ve_k+ve;
    
    u(i)= -Ka(1:4)*xobs - Ka(5)*ve_k; % accion de control con observador
    
%     u(i)= -K*xobs + G*ref; % accion de control con observador
    
    % No linealidad del controlador
    alin = 0.001;
    if (abs(u(i))<alin)
        acc(i) = 0;
    elseif (u(i)>alin)
        acc(i) = u(i)-alin;
    elseif (u(i)<-alin)
        acc(i) = u(i)+alin;
    end
    
    if (acc(i)>1)
        acc(i) = 1;
    end
    
    if (acc(i)<-1)
        acc(i) = -1;
    end

%     x = modelo_avion(Ts,x,u(i),config); % simulacion del comportamiento real del avion
    x = modelo_avion(Ts,x,acc(i),config); % simulacion del comportamiento real del avion
    
%     xobs = A*xobs + B*u(i) + Ko'*(Y-Yobs); % nuevo estado del observador
    xobs = A*xobs + B*acc(i) + Ko'*(Y-Yobs); % nuevo estado del observador
    
    Yobs = C*xobs; % salida con estado del observador
    
    alpha(i)=x(1);
    phi(i)=x(2);
    phi_p(i)=x(3);
    h(i)=x(4);
    
    alpha_obs(i)=xobs(1);
    phi_obs(i)=xobs(2);
    phi_p_obs(i)=xobs(3);
    h_obs(i)=xobs(4);
    
end

%% Plots

% Sin observador

figure(1);

subplot(4,2,1);
plot(t_v,alpha);
grid on; 
title('\alpha');
hold on;

subplot(4,2,2);
plot(t_v,phi);
grid on;
title('\phi');
hold on;

subplot(4,2,3); 
plot(t_v,phi_p);
grid on;
title('\phi_p');
hold on;

subplot(4,2,4);
plot(t_v,h);
grid on;
title('Altura h');
hold on;

subplot(4,1,3);
plot(t_v,u);
grid on;
title('Acción de control');
xlabel('Tiempo [s]');

subplot(4,1,4);
plot(t_v,acc);
grid on;
title('Acción de control (acc)');
xlabel('Tiempo [s]');

% con Observador

% figure(2);
% 
% subplot(4,2,1);
% plot(t_v,alpha_obs);
% grid on; 
% title('\alpha');
% hold on;
% 
% subplot(4,2,2);
% plot(t_v,phi_obs);
% grid on;
% title('\phi');
% hold on;
% 
% subplot(4,2,3); 
% plot(t_v,phi_p_obs);
% grid on;
% title('\phi_p');
% hold on;
% 
% subplot(4,2,4);
% plot(t_v,h_obs);
% grid on;
% title('Altura h');
% hold on;
% 
% subplot(4,1,3);
% plot(t_v,u);
% grid on;
% title('Acción de control');
% % xlabel('Tiempo [s]');
% 
% subplot(4,1,4);
% plot(t_v,acc);
% grid on;
% title('Acción de control limitada');
% xlabel('Tiempo [s]');


% Nota: controlador poco robusto ante no linealidades
% probar con un controlador con integrador

