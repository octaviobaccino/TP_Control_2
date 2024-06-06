clear; close all; clc

%--------------------------------------------------------------------------
% Alumno: Octavio Baccino
% Profesor: Julian Pucheta
% Trabajo Práctico 2
%--------------------------------------------------------------------------
% Caso de estudio 3. Sistema no lineal de cuatro variables de estado
%--------------------------------------------------------------------------
% Ítem [3] 
% Calcular sistema controlador que haga evolucionar al péndulo en el equilibrio estable. 
% Objetivo de control: partiendo de una condición inicial nula en el desplazamiento y el ángulo en pi, hacer 
% que el carro se desplace a 10 metros evitando las oscilaciones de la masa m, considerando que es una 
% grúa. Una vez que ?=10 modificar a m a un valor 10 veces mayor y volver al origen evitando oscilaciones.

%% Matrices del sistema en variables de estado

% x = [delta; delta_p; phi; phi_p];

m = 0.1; 
F = 0.1; 
L = 1.6; 
g = 9.8; 
M = 1.5;

T0 = 1e-4; % tiempo integracion de euler

% Matrices del sistema continuo

A = [0 1 0 0;
    0 -F/M -m*g/M 0;
    0 0 0 1;
    0 -F/L/M -g*(m+M)/L/M 0];

B = [0; 1/M; 0; 1/L/M];

C = [1 0 0 0
    0 0 1 0];

D = 0;

% Matrices ampliadas

Aa = [A zeros(4,1);
    -C(1,:) 0];

Ba = [B; 0];

% Calculo del regulador lineal cuadratico (LQR)

Q = diag([1e-1 1e1 1e1 1e3 1e2/2]);

R = 100;

[Ka,S,P] = lqr(Aa,Ba,Q,R);

fprintf('Matriz K ampliada:');
Ka

%% Observador sistema inicial

Ao = A';
Bo = C';
Co = B';
Do = D';

% Calculo del controlador para el observador

Qo = diag([1e1 1e2 1e5 1e6]);

Ro = diag([1e-0 1e-0]);

[Ko,So,Po] = lqr(Ao,Bo,Qo,Ro);

fprintf('Matriz Ko del observador:');
Ko

%% Sistema con cambio de masa

m1 = 10*m; 


% Matrices del sistema continuo

A1 = [0 1 0 0;
    0 -F/M -m1*g/M 0;
    0 0 0 1;
    0 -F/L/M -g*(m1+M)/L/M 0];

B1 = [0; 1/M; 0; 1/L/M];

C1 = [1 0 0 0
    0 0 1 0];

D1 = 0;

% Matrices ampliadas

A1a = [A1 zeros(4,1);
    -C1(1,:) 0];

B1a = [B; 0];

% Calculo del regulador lineal cuadratico (LQR)

Q1 = diag([1e-1 1e1 1e1 1e3 1e2/2]);

R1 = 100;

[Ka1,S1,P1] = lqr(A1a,B1a,Q1,R1);

fprintf('Matriz K ampliada:');
Ka1

%% Observador sistema con cambio de masa

A1o = A1';
B1o = C1';
C1o = B1';
D1o = D1';

% Calculo del controlador para el observador

Q1o = diag([1e1 1e5 1e3 1e7]);

R1o = diag([1e-0 1e-0]);

[K1o,S1o,P1o] = lqr(A1o,B1o,Q1o,R1o);

fprintf('Matriz Ko del observador:');
K1o


%% Simulacion

tsim = 100;
t_v = 0:T0:tsim-T0;

delta_v = zeros(1,length(t_v));
delta_p_v = zeros(1,length(t_v));
phi_v = zeros(1,length(t_v));
phi_p_v = zeros(1,length(t_v));

delta_o_v = zeros(1,length(t_v));
delta_p_o_v = zeros(1,length(t_v));
phi_o_v = zeros(1,length(t_v));
phi_p_o_v = zeros(1,length(t_v));

u = zeros(1,length(t_v));
acc = zeros(1,length(t_v));

psi = zeros(1,length(t_v));

y = zeros(1,length(t_v));

ref = [10*ones(1,length(t_v)/2) zeros(1,length(t_v)/2)];

x = [0; 0; pi; 0]; % condiciones iniciales

xop=[0; 0; pi; 0]; % punto de operacion

integracion = 0;

delta_pp = 0; % condiciones iniciales
phi_pp = 0;

delta = 0;
delta_p = 0;
phi = pi;
phi_p = 0;

Yobs = [0; pi];

xobs = [0; 0; pi; 0];

% estados = [0;0;pi;0];

for i = 1:length(t_v)
    
    Y = C*(x);
    y(i) = Y(1,1);
    
    psi_p = ref(i) - y(i);
    psi(i) = integracion + psi_p*T0;
    
    if (ref(i) == 10)
        u(i)=-Ka(1:4)*(xobs-xop)-Ka(5)*psi(i); % agrego xop para que la accion inicial no sea tan grande
%         u(i)=-Ka(1:4)*(x-xop)-Ka(5)*psi(i);
    end
    
    if (ref(i) == 0)
        u(i)=-Ka1(1:4)*(xobs-xop)-Ka1(5)*psi(i); % agrego xop para que la accion inicial no sea tan grande
%         u(i)=-Ka1(1:4)*(x-xop)-Ka1(5)*psi(i);
    end
    
    % Alienalidad del controlador
    alin = 0.15;
    if (abs(u(i))<alin)
        acc(i) = 0;
    elseif (u(i)>alin)
        acc(i) = u(i)-alin;
    elseif (u(i)<-alin)
        acc(i) = u(i)+alin;
    end
    
    % Modelo no lineal (real)
    delta = x(1);
    delta_p = x(2);
    phi = x(3);
    phi_p = x(4);
    
    if (ref(i) == 10)
        % controlador normal
%         delta_pp = (1/(M+m))*(u(i)-m*L*phi_pp*cos(phi)+m*L*phi_p^2*sin(phi)-F*delta_p);
%         phi_pp = (1/L)*(g*sin(phi)-delta_pp*cos(phi));
        % controlador con alinealidad
        delta_pp = (1/(M+m))*(acc(i)-m*L*phi_pp*cos(phi)+m*L*phi_p^2*sin(phi)-F*delta_p);
        phi_pp = (1/L)*(g*sin(phi)-delta_pp*cos(phi));
    end
    
    if (ref(i) == 0)
        % controlador normal
%         delta_pp = (1/(M+m1))*(u(i)-m1*L*phi_pp*cos(phi)+m1*L*phi_p^2*sin(phi)-F*delta_p);
%         phi_pp = (1/L)*(g*sin(phi)-delta_pp*cos(phi));
        % controlador con alinealidad
        delta_pp = (1/(M+m1))*(acc(i)-m1*L*phi_pp*cos(phi)+m1*L*phi_p^2*sin(phi)-F*delta_p);
        phi_pp = (1/L)*(g*sin(phi)-delta_pp*cos(phi));
    end
    
    xp = [delta_p; delta_pp; phi_p; phi_pp];
    
    x = x + T0*xp;
    
    delta_v(i) = x(1);
    delta_p_v(i) = x(2);
    phi_v(i) = x(3);
    phi_p_v(i) = x(4);
    
    % Observador-----------------------------------------------------------
    
%     if (ref(i) == 10)
%         xobs_p = A*(xobs-xop) + B*u(i) + Ko'*(Y-Yobs);
%         xobs = xobs + xobs_p*T0;
%     end
%     
%     if (ref(i) == 0)
%         xobs_p = A*(xobs-xop) + B*u(i) + K1o'*(Y-Yobs);
%         xobs = xobs + xobs_p*T0;
%     end
    
    % Para no linealidad del controlador
    if (ref(i) == 10)
        xobs_p = A*(xobs-xop) + B*acc(i) + Ko'*(Y-Yobs);
        xobs = xobs + xobs_p*T0;
    end
    
    if (ref(i) == 0)
        xobs_p = A*(xobs-xop) + B*acc(i) + K1o'*(Y-Yobs);
        xobs = xobs + xobs_p*T0;
    end
    
    delta_o_v(i) = xobs(1);
    delta_p_o_v(i) = xobs(2);
    phi_o_v(i) = xobs(3);
    phi_p_o_v(i) = xobs(4);
    
    Yobs = C*(xobs);
    
    % Modelo lineal
%     xp = A*(x-xop) + B*u(i);
%     x = x + xp.*T0;
%     
%     delta_v(i) = x(1);
%     delta_p_v(i) = x(2);
%     phi_v(i) = x(3);
%     phi_p_v(i) = x(4);
    
    integracion=psi(i);
    
end

%% Plots

% Sin observador

figure(1);

subplot(3,2,1);
plot(t_v,delta_v)
hold on;
plot(t_v,ref,'k--')
plot(t_v,delta_o_v,'--')
grid on; 
title('\delta');
hold on;

subplot(3,2,2);
plot(t_v,delta_p_v);
hold on;
plot(t_v,delta_p_o_v,'--');
grid on;
title('\delta_p');
hold on;

subplot(3,2,3); 
plot(t_v,phi_v);
hold on;
plot(t_v,phi_o_v,'--');
grid on;
title('\phi');
hold on;

subplot(3,2,4);
plot(t_v,phi_p_v);
hold on;
plot(t_v,phi_p_o_v,'--');
grid on;
title('\phi_p');
legend('sin Observador','con Observador')
hold on;

subplot(3,1,3);
plot(t_v,u);
grid on;
title('Accion de control');
hold on;

% Preguntas: 

% No linealidades hacen oscilar el sistema?

% Como ajustar controlador del avion



