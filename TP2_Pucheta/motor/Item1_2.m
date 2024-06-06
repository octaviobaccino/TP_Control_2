clear; close all; clc

%--------------------------------------------------------------------------
% Alumno: Octavio Baccino
% Profesor: Julian Pucheta
% Trabajo Práctico 2
%--------------------------------------------------------------------------
% Caso de estudio 1. Sistema de tres variables de estado
%--------------------------------------------------------------------------
% -Asumiendo que no puede medirse directamente la corriente, pero sí la 
% velocidad y el ángulo, proponer un controlador que logre el objetivo. 
% -Determinar el efecto de la nolinealidad en la acción de control, 
% descripta en la Fig. 2, y verificar cuál es el máximo valor admisible de ésa no linealidad. 

%% Disenio de controlador con observador

% Importar valores de Excel

[t_v, w_v, i_v, Va_v, TL_v] = importfile_MOTOR("Curvas_Medidas_Motor_2024.xls"); 

% Curvas de referencia del motor
if false
    figure;

    subplot(4,1,1)
    plot(t_v,w_v)

    subplot(4,1,2)
    plot(t_v,i_v)

    subplot(4,1,3)
    plot(t_v,Va_v)

    subplot(4,1,4)
    plot(t_v,TL_v)
end

% Modelado del motor en variables de estado
Ra = 28.13;
Laa = 5.7586e-04;
Km = 0.0605;
Ki = 0.0122;
J = 2e-9;
Bm = 0;

% Matrices del sistema
A = [-Ra/Laa -Km/Laa 0 ; Ki/J -Bm/J 0 ; 0 1 0];
B = [1/Laa 0 0]';
C = [0 1 0; 0 0 1]; % 2 Salidas y = [w;theta]
% D = 0;

% Matrices ampliadas
Aa = [A zeros(3,1); -C(2,:) 0];
Ba = [B; 0;];
Ca = [C zeros(2,2)];

% LQR
Q=1*diag([1e1 1e-2 1e0 1e6]);    R=1e-1;
Kamp=lqr(Aa,Ba,Q,R);

fprintf('Polos de lazo cerrado:');
eig(Aa-Ba*Kamp)

%% Observador
% Propongo un controlador con obervador ya que no se puede medir la
% corriente

Ao = A';
Bo = C';
Co = B';
D = 0;

Qo=diag([1e1 1e0 1e2]);    Ro=diag([1e-5 1e-6]);

Ko=lqr(Ao,Bo,Qo,Ro);

% Vectores de tiempo
tsim=21; 
Ts=1e-5; 
% t_v=0:Ts:(tsim-Ts);

delay = 1;

t_delay = 0:Ts:delay;
t_1 = t_delay(end)+Ts:Ts:6;
t_2 = t_1(end)+Ts:Ts:11;
t_3 = t_2(end)+Ts:Ts:16;
t_4 = t_3(end)+Ts:Ts:tsim;

t_v = [t_delay t_1 t_2 t_3 t_4];

% Vector de referenia
ref = [zeros(1,length(t_delay)) pi/2*ones(1,length(t_1)) ...
    -pi/2*ones(1,length(t_2)) pi/2*ones(1,length(t_3)) -pi/2*ones(1,length(t_4))];

% Vector de torque
TL_sim = [zeros(1,length(t_delay)) zeros(1,length(t_1)/2) 1e-3+1e-5*randn(1,length(t_1)/2) ...
    zeros(1,length(t_2)/2) 1e-3+1e-5*randn(1,length(t_2)/2) ...
    zeros(1,length(t_3)/2) 1e-3+1e-5*randn(1,length(t_3)/2) ...
    zeros(1,length(t_4)/2) 1e-3+1e-5*randn(1,length(t_4)/2)];

% TL_sim = [zeros(1,length(t_delay)) zeros(1,length(t_1)/2) 1e-3+1e-5*ones(1,length(t_1)/2) ...
%     zeros(1,length(t_2)/2) 1e-3+1e-5*ones(1,length(t_2)/2) ...
%     zeros(1,length(t_3)/2) 1e-3+1e-5*ones(1,length(t_3)/2) ...
%     zeros(1,length(t_4)/2) 1e-3+1e-5*ones(1,length(t_4)/2)];

%condiciones iniciales

ia = zeros(1,length(t_v));
omega = zeros(1,length(t_v));
theta = zeros(1,length(t_v));

ia(1)=0;                %Corriente de armadura          
omega(1)=0;             %Valocidad angular 
theta(1)=0;             %Posicion angular 

estados=[ia(1);
        omega(1);
        theta(1)];
    
ia_o = zeros(1,length(t_v));
omega_o = zeros(1,length(t_v));
theta_o = zeros(1,length(t_v));

ia_o(1)=0;                %Corriente de armadura          
omega_o(1)=0;             %Valocidad angular 
theta_o(1)=0;             %Posicion angular 
    
estados_obs=[ia_o(1);
            omega_o(1);
            theta_o(1)];
    
% Xop=[0 0 0]';
x=[ia(1) omega(1) theta(1)]'; % vector de estados
xobs=[0 0 0]'; %inicializacion para el observador

acc = zeros(1,length(t_v));
u = zeros(1,length(t_v));
psi = zeros(1,length(t_v));

psi(1)=0;
integracion(1)=psi(1);

y = zeros(1,length(t_v));
y_o = zeros(1,length(t_v));

for i=1:length(t_v)
    
    %------sin Observador----------------------
    
    Y = C*x;
    y(i)=Y(2,1);
    
    psi_p=ref(i) - y(i);
    psi(i)=integracion+psi_p*Ts;
    
    u(i) = -Kamp(1:3)*xobs-Kamp(4)*psi(i);
%     u(i) = -Kamp(1:3)*xobs-Kamp(4)*psi(i); % u con x del observador

    % No linealidad del controlador
    alin = 1;
    if (u(i)<alin) && (u(i)>-alin)
        acc(i) = 0;
    elseif (u(i)>alin)
        acc(i) = u(i)-alin;
    elseif (u(i)<-alin)
        acc(i) = u(i)+alin;
    end
    
    %Variables del sistema lineal
    ia(i)= x(1);
    omega(i)= x(2);
    theta(i)= x(3);
    
%     x1_p=-Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x1_p=-Ra*x(1)/Laa-Km*x(2)/Laa+acc(i)/Laa;
    x2_p=Ki*x(1)/J-Bm*x(2)/J-TL_sim(i)/J;
    x3_p=x(2);
    
    xp=[x1_p; x2_p; x3_p];
    x=x+Ts*xp;
    
%     estados=[ia(i);omega(i);theta(i)];
%     integracion=psi(i);
    
    %------con Observador----------------------

    Y_o = C*xobs;
    y_o(i) = Y_o(2,1);
    
%     xobs_p = A*xobs + B*u(i) + Ko'*(Y-Y_o);
    xobs_p = A*xobs + B*acc(i) + Ko'*(Y-Y_o);
    
    xobs = xobs + xobs_p*Ts;
    
    ia_o(i)= xobs(1);
    omega_o(i)= xobs(2);
    theta_o(i)= xobs(3);
    
%     %--------------------------------------------
    
%     estados=[ia(i);omega(i);theta(i)];
    integracion=psi(i);
%     estados_obs=[ia_o(i);omega_o(i);theta_o(i)];
end

%% Plots

% Plots control con integracion
figure(2)

subplot(3, 2, 1);
hold on
plot(t_v,ia);
hold off
title('Corriente de armadura i_a');
xlabel('Tiempo [s]');
ylabel('Corriente [A]');
grid on;

subplot(3, 2, 2);
hold on
plot(t_v,omega);
hold off
title('Velocidad angular \omega_r');
xlabel('Tiempo [s]');
ylabel('Velocidad angular [rad/s]');
grid on;

subplot(3, 2, 3);
hold all
plot(t_v,theta);
plot(t_v,ref,'--');
hold off
title('Poscion angular \theta_t');
xlabel('Tiempo [s]');
ylabel('Posicion angular [Rad]');
legend('\theta','ref');
grid on;

subplot(3, 2, 4);
hold all
plot(t_v,u);
% plot(t_v,acc);
hold off
title('Accion de control u_t');
xlabel('Tiempo [s]');
ylabel('V');
grid on;

subplot(3, 1, 3);
hold all
plot(t_v,acc);
% plot(t_v,acc);
hold off
title('Accion de control acc');
xlabel('Tiempo [s]');
ylabel('V');
grid on;

figure(3)

subplot(2, 1, 1);
hold on
plot(t_v,TL_sim);
hold off
title('Torque de entrada TL');
xlabel('Tiempo [s]');
ylabel('Nm');
grid on;

subplot(2, 1, 2);
hold on
plot(omega(1:round(length([t_delay t_1])/2)),theta(1:round(length([t_delay t_1])/2)));
hold off
title('Vel. Angular vs Angulo');
xlabel('rad/s');
ylabel('rad');
grid on;

% Plots con observador
figure(4)

subplot(2,1,1);
hold all
plot(t_v,theta);
plot(t_v,theta_o);
hold off
title('Poscion angular \theta_t');
xlabel('Tiempo [s]');
ylabel('\theta_t [rad]');
legend('Sin obervador','Con observador');
grid on;

subplot(2,1,2);
hold all
plot(t_v,omega);
plot(t_v,omega_o);
hold off
title('Velocidad angular \omega_t');
xlabel('Tiempo [s]');
ylabel('\omega_t [rad/s]');
legend('Sin obervador','Con observador');
grid on;


