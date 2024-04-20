clear; close all;

%--------------------------------------------
% Ejercicio 5 Trabajo Practico Pucheta
%--------------------------------------------

%% Metodo de la respuesta al escalon

% Importar valores de Excel
[t_v, w_v, i_v, Va_v, TL_v] = importfile_MOTOR("Curvas_Medidas_Motor_2024.xls"); 

% plots para visualizar comprotamiento del motor

enable_plot = 0;

if enable_plot == 1
    
    figure
    subplot(4,1,1)
    plot(t_v, Va_v)
    subplot(4,1,2)
    plot(t_v, i_v)
    subplot(4,1,3)
    plot(t_v, w_v)
    subplot(4,1,4)
    plot(t_v, TL_v)
    
end

%% Obtencion de la funcion de transferencia aproxiamda con metodo de Chen

% Valores de t1,t2,t3 | y1,y2,y3 | delay
% obtenidos observando el grafico

delay = 0.0351;

t_sel = 0.000075;

t1 = t_sel;
t2 = 2*t_sel;
t3 = 3*t_sel;

y1 = interp1(t_v,w_v,t1+delay);
y2 = interp1(t_v,w_v,t2+delay);
y3 = interp1(t_v,w_v,t3+delay);

StepAmplitude = 12;

K = 198.2/12;

k1=(1/StepAmplitude)*y1/K-1; %Afecto el valor del Escalon

k2=(1/StepAmplitude)*y2/K-1;

k3=(1/StepAmplitude)*y3/K-1;

be=4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;

if be > 0
    
    alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
    alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));

else
    
    alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
    alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));

end

beta=(k1+alfa2)/(alfa1-alfa2);

T1 = -t1/log(alfa1);
T2 = -t1/log(alfa2);
T3 = beta*(T1-T2)+T1;

G_ang = tf(K*[T3 1],conv([T1 1],[T2 1]));

% step(G_ang) --> para verificar dinamica del sistema

% G_ang =
%  
%        -0.0001741 s + 16.52
%   -------------------------------
%   1.563e-09 s^2 + 7.635e-05 s + 1

% ignoro el 0 de la TF

% G = 1/(1.563e-09*s^2 + 7.635e-05*s + 1)

% eqn1 =  Ki/(Bm*Ra+Km*Ki) == 16.52

% eqn2 = J*Laa/(Bm*Ra+Km*Ki) == 1.563e-09

% eqn3 = (Bm*Laa+Ra*J)/(Bm*Ra+Km*Ki) == 7.635e-05

% Ra = Va/Imax = 12/0.4266 = 28.13

% Bm = 0 ==> Ki/(Ki*Km) = 16.52 ==> 1/Km = 16.52 ==> Km = 0.0605

% Bm = 0 ==> (Ra*J)/(Km*Ki) == 7.635e-05

% Me quedan dos ecuaciones y tres incognitas, por lo tanto elijo un valor
% coherente de alguna de las varialbes y voy cambiandolo manualmente para
% ajustar su valor a uno adecuado

% Laa = 1.563e-09/7.635e-05*Ra

% Ki = J*Laa/Km/1.563e-09

Ra = 28.13;
Km = 0.0605;
% Laa = 366e-6;
% Laa = 1.563e-09/7.635e-05*Ra;
Laa = 5.7586e-04;
J = 2e-9; % este valor fue elgido por mi
Ki = J*Laa/Km/1.563e-09;
B = 0;

%% Evolucion del sistema

A_m = [-Ra/Laa -Km/Laa 0 ; Ki/J -B/J 0 ; 0 1 0];

B1_m = [1/Laa 0 0]'; % reemplazar 1 por Va y TL 

B2_m = [0 -1/J 0]';

x = [0 0 0]';

Ts = 1e-5;

t_delay = 0:Ts:delay;
 
t_1 = t_delay(end)+Ts:Ts:0.1869;

t_2 = t_1(end)+Ts:Ts:0.3366;

t_3 = t_2(end)+Ts:Ts:0.4872;

t_4 = t_3(end)+Ts:Ts:0.6;

t_line = [t_delay t_1 t_2 t_3 t_4];

Va_sim = [zeros(1,length(t_delay)) 12*ones(1,length([t_1 t_2 t_3 t_4]))];

TL_sim = [zeros(1,length([t_delay t_1])) 1e-3+1e-5*randn(1,length(t_2)) ...
          zeros(1,length(t_3)) 1e-3+1e-5*randn(1,length(t_4))];

% u = [Va_sim; TL_sim; zeros(1,length(t_line))]';      
      
ia = zeros(length(t_v),1);
wr = zeros(length(t_v),1);
theta = zeros(length(t_v),1);

for idx = 1:(length(t_line))
    
    ia(idx) = x(1);
    wr(idx) = x(2);
    theta(idx) = x(3);
    
    xp = A_m*x + B1_m.*Va_sim(idx) + B2_m.*TL_sim(idx); % agragar "u" para variar las condiciones
    x = x + xp.*Ts;
    
end

%% Plots para comparar

fz = 11;

figure

subplot(4,1,1)
plot(t_line, Va_sim, 'LineWidth', 1.2)
hold on;
plot(t_v, Va_v,'--', 'LineWidth', 1.2)
title('Tension de entrada (Va)','Interpreter','latex','FontSize',fz)
% xlabel('t[s]','Interpreter','latex','FontSize',fz)
ylabel('V[V]','Interpreter','latex','FontSize',fz)

subplot(4,1,2)
plot(t_line, ia, 'LineWidth', 1.2)
hold on;
plot(t_v, i_v, '--', 'LineWidth', 1.2)
title('Corriente (Ia)','Interpreter','latex','FontSize',fz)
% xlabel('t[s]','Interpreter','latex','FontSize',fz)
ylabel('I[A]','Interpreter','latex','FontSize',fz)

subplot(4,1,3)
plot(t_line, wr, 'LineWidth', 1.2)
hold on;
plot(t_v, w_v,'--', 'LineWidth', 1.2)
title('Velocidad angular (Wr)','Interpreter','latex','FontSize',fz)
% xlabel('t[s]','Interpreter','latex','FontSize',fz)
ylabel('w[rad/s]','Interpreter','latex','FontSize',fz)

subplot(4,1,4)
plot(t_line,TL_sim, 'LineWidth', 1.2)
hold on;
plot(t_v, TL_v,'--', 'LineWidth', 1.2)
title('Perturbaciones (TL)','Interpreter','latex','FontSize',fz)
xlabel('t[s]','Interpreter','latex','FontSize',fz)
ylabel('TL[Nm]','Interpreter','latex','FontSize',fz)
legend('Identificada','Medida','Interpreter','latex','FontSize',fz-1)

