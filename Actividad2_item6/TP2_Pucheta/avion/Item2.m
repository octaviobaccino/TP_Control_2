clear; close all; clc

%--------------------------------------------------------------------------
% Alumno: Octavio Baccino
% Profesor: Julian Pucheta
% Trabajo Práctico 2
%--------------------------------------------------------------------------
% Caso de estudio 2. Sistema de cuatro variables de estado
%--------------------------------------------------------------------------
% Ítem [2]
% Para el caso del avión, emplear un tiempo de integración por Euler 
% adecuado y un tiempo de simulación de 70seg. Los parámetros son a=0.07; 
% ?=9; b=5; c=150, hallar un controlador para que los polos de lazo cerrado 
% se ubican en ?i=-15?15j; -0.5?0.5j, para referencias de 100 y -100 
% metros en altura, ambas con alturas iniciales de -500 y 500. 
% -Proponer un controlador en tiempo discreto en variables de estado para 
% que el proceso evolucione en los rangos de validez del modelo, es decir 
% donde los ángulos y el valor de la acción de control en valor absoluto 
% son menores a la unidad. 
%--------------------------------------------------------------------------

%% Matrices para el sistema en variables de estado

% Controlador por asigancion de polos

% x = [alpha; phi; phi_p; h]

a = 0.07;
b = 5;
w = 9;
c = 150;

A = [-a a 0 0; % Matriz de estados
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];

B = [0; 0; w^2*b; 0]; % Matriz de entradas

C = [0 0 0 1]; % Matriz de salidas

D = 0; % Matriz de transmision directa

M=[B A*B A^2*B A^3*B]; %Matriz Controlabilidad

%Cálculo del controlador por asignación de polos
auto_val = eig(A);
c_ai=poly(auto_val);

W=[c_ai(4) c_ai(3) c_ai(2) 1;
 c_ai(3) c_ai(2) 1 0;
 c_ai(2) 1 0 0;
 1 0 0 0];

T = M*W;

A_controlable=inv(T)*A*T; %Verificación de que T esté bien

%Ubicación de los polos de lazo cerrado en mui:
mui(1) = -15 + 15j; 
mui(2) = conj(mui(1)); 
mui(3) = -0.5 + 0.5i;
mui(4) = conj(mui(3));

alfa_ia = poly(mui);

K=fliplr((alfa_ia(2:5)-c_ai(2:5)))*inv(T);

Gj=-inv(C*inv(A-B*K)*B); % Ganancia para llegar a la referencia
eig(A-B*K) 
%Verifico que todos los polos estén en el semiplano izquierdo
% polo mas alejado del origen (dominante) -21, por lo tanto el tiempo 
% asociado es log(0.95)/21 = 2e-3
Ts = 2e-4;

delay = 10;

t_delay = 0:Ts:delay;
t1 = t_delay(end):Ts:45;
t2 = t1(end):Ts:80;

t_v = [t_delay t1 t2];

ref = [-500*ones(1,length(t_delay)) 100*ones(1,length(t1)) -100*ones(1,length(t2))];

x = [0; 0; 0; -500]; % condiciones iniciales

alpha = zeros(1,length(t_v));
phi = zeros(1,length(t_v));
phi_p = zeros(1,length(t_v));
h = zeros(1,length(t_v));
u = zeros(1,length(t_v));
y = zeros(1,length(t_v));

for i = 1:length(t_v)
    
    u(i)=-K*x+Gj*ref(i);
    
    xp = A*x + B*u(i);
    x = x + xp.*Ts;
    
    y(i) = C*x;
    
    alpha(i) = x(1);
    phi(i) = x(2);
    phi_p(i) = x(3);
    h(i) = x(4);
    
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

% figure(2);hold on;
% subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
% subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;

