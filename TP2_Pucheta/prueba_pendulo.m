clc; clear all; close all;

% Ítem [6] Incorporar un observador para el caso en que sólo puedan medirse el desplazamiento delta
% y el ángulo theta, repetir las simulaciones para las condiciones anteriores y graficar los resultados en 
% gráficas superpuestas para el equilibrio estable. Con punto 3 bien hecho


m=.1;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;

% Matrices del sistema linealizado equilibrio estable
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
B=[0; 1/M; 0; 1/(l*M)];
C=[1 0 0 0]; 
%usar modelo linealizado 180/pi


% Construcción del sistema ampliado
Aa=[A zeros(4,1); -C 0];
Ba=[B; 0];
Ca=[C 0];

%Diseño con LQR
% Q=1*diag([10 .01 10 1]);    R=5;
Q=1*diag([0.1 100 10 10000 1]);    R=100;
Kamp=lqr(Aa,Ba,Q,R); %controlador para primera etapa hasta 2 metros sin carga
%ampliado por el integrador
K1=Kamp(1:4);
K1I=-Kamp(5);


%Sistema con 10 veces mas masa---------------------------------------------
m=.1*10;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;
Am=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
Bm=[0; 1/M; 0; 1/(l*M)];
Cm=[1 0 0 0]; 

% Construcción del sistema ampliado
Am=[Am zeros(4,1); -Cm 0];
Bm=[Bm; 0];
Cm=[C 0];

%Diseño con LQR
% Q=1*diag([10 .01 10 1]);    R=5;
Qm=1*diag([1 1000 1 10000 1]);    Rm=100;
Kamp_m=lqr(Am,Bm,Qm,Rm); %controlador para primera etapa hasta 2 metros sin carga
%ampliado por el integrador
K1m=Kamp_m(1:4);
K1Im=-Kamp_m(5);
%----------------------Observador------------------------------------------
C=[1 0 0 0;0 0 1 0];
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];

%primer masa
mo1=0.1*1;
A1=[0 1 0 0;0 -Fricc/M -mo1*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(mo1+M)/(l*M) 0];

Ao1=A1';
Bo1=C';

Qo1=1*diag([1 1/100 1 1]);    Ro1=diag([2 1]); 
Ko1=lqr(Ao1,Bo1,Qo1,Ro1);

Ko1=place(Ao1,Bo1,[-101 -200 -300 -400])


%carga de masa segunda 1kg
mo2=0.1*10;
A2=[0 1 0 0;0 -Fricc/M -mo2*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(mo2+M)/(l*M) 0];

Ao2=A2';
Bo2=C';

Qo2=1*diag([1 1/10 1 1]);    Ro2=diag([1 1]); 
Ko2=lqr(Ao2,Bo2,Qo2,Ro2);


%--------------------------------------------------------------------------

%Simulación del control:
h=10e-5;%paso
tsim=100; %tiempo de simulacion
t=0:h:(tsim-h);

%Referencia
% sp_dist=2; % la distancia de desplazamiento es -10m 
%Nueva referencia que en 10 segundos va a 2m y despues vuelve a 0m en 10s
sp_dist=1*square(2*pi*t/tsim)+1;
ref_ang=pi;
pasos=round(tsim/h);
m=1;
m=ones(1,pasos);
m=m*0.1;
m((pasos/2):end)=m((pasos/2):end)*10;

% figure(2)
% plot(t,m)

%condiciones iniciales
delta(1)=0;        %x1
delta_p(1)=0;      %x2
theta(1)=pi;     %x3
theta_p(1)=0;      %x4
psi(1)=0;

estados=[delta(1);
        delta_p(1);
        theta(1);
        theta_p(1)];
integracion(1)=psi(1);    
Xop=[0 0 pi 0]';
x=[delta(1) delta_p(1) theta(1) theta_p(1)]';

theta_pp=0;

%inicializo los observadores
estados_obs=[delta(1);
    delta_p(1);
    theta(1);
    theta_p(1)];

xobs=[delta(1) delta_p(1) theta(1) theta_p(1)]';

%Inicializo los controladores
K=K1;KI=K1I;Ko=Ko1;
for i=1:round(tsim/h)
    
    if m(i)>=0.5
        K=K1m;
        KI=K1Im;
        Ko=Ko2;
    end
    
    psi_p=sp_dist(i)-C(1,:)*(estados_obs-Xop);
    psi(i)=integracion+psi_p*h;
       
    u(i)=-K*(estados_obs-Xop)+psi(i)*KI; 
    
    %Variables del sistema lineal
    delta(i)= x(1);
    delta_p(i)= x(2);
    theta(i)= x(3);
    theta_p(i)= x(4);
    
    
    delta_pp=(u(i)-Fricc*x(2)-m(i)*l*theta_pp*cos(x(3))+m(i)*l*sin(x(3))*x(4)^2)/(M+m(i));
    theta_pp=((g*sin(x(3)))-delta_pp*cos(x(3)))/l;%no restarle 
    x_p_1=x(2);
    x_p_2=delta_pp;
    x_p_3=x(4);
    x_p_4=theta_pp;
    xp=[x_p_1;x_p_2;x_p_3;x_p_4];
    x=x+h*xp;
    %Observador
    delta_o(i)= xobs(1);
    delta_p_o(i)= xobs(2);
    theta_o(i)= xobs(3);
    theta_p_o(i)= xobs(4);
    
    
    y_sal_o = C * (estados_obs);%2x4 * 4*1=2x1
    y_sal   = C * (estados);%2x1
    
    e=y_sal-y_sal_o;
    
    x_antp     = A*(xobs-Xop)+B*u(i)+Ko'*e;%complicacion
    xobs       = xobs + x_antp*h;
    
    %----------------------------
    estados=[delta(i);
        delta_p(i);
        theta(i);
        theta_p(i)];
    
        estados_obs=[delta_o(i);
        delta_p_o(i);
        theta_o(i);
        theta_p_o(i)];
    
    integracion=psi(i);
    
end
    
    
%_______________SIN OBS---------------------------------------------------

%condiciones iniciales
delta_so(1)=0;        %x1
delta_p_so(1)=0;      %x2
theta_so(1)=pi;     %x3
theta_p_so(1)=0;      %x4
psi_so(1)=0;

estados_so=[delta_so(1);
        delta_p_so(1);
        theta_so(1);
        theta_p_so(1)];
integracion_so(1)=psi_so(1);    
Xop_so=[0 0 pi 0]';
x_so=[delta_so(1) delta_p_so(1) theta_so(1) theta_p_so(1)]';

theta_pp_so=0;

%Inicializo los controladores
K=K1;KI=K1I;
for i=1:round(tsim/h)
    
    if m(i)>=0.5
        K=K1m;
        KI=K1Im;
    end
    
    psi_p_so=sp_dist(i)-C(1,:)*(estados_so-Xop_so);
    psi_so(i)=integracion_so+psi_p_so*h;
    
    u_so(i)= -K*(estados_so(1:4)-Xop_so(1:4))+KI*psi_so(i);
    
    %Variables del sistema lineal
    delta_so(i)= x_so(1);
    delta_p_so(i)= x_so(2);
    theta_so(i)= x_so(3);
    theta_p_so(i)= x_so(4);
    
    
    delta_pp_so=(u_so(i)-Fricc*x_so(2)-m(i)*l*theta_pp_so*cos(x_so(3))+m(i)*l*sin(x_so(3))*x_so(4)^2)/(M+m(i));
    theta_pp_so=((g*sin(x_so(3)))-delta_pp_so*cos(x_so(3)))/l;%no restarle 
    x_p_1_so=x_so(2);
    x_p_2_so=delta_pp_so;
    x_p_3_so=x_so(4);
    x_p_4_so=theta_pp_so;
    xp_so=[x_p_1_so;x_p_2_so;x_p_3_so;x_p_4_so];
    x_so=x_so+h*xp_so;
    
    estados_so=[delta_so(i);
        delta_p_so(i);
        theta_so(i);
        theta_p_so(i)];
    integracion_so=psi_so(i);
    
end
    
%----------------------------------------------------------------------


figure
subplot(3, 2, 1);
hold on
plot(t,delta);
plot(t,delta_so);
plot(t,sp_dist,'g--');
hold off
legend({'Con observador','Sin observador','Referencia'})
title('desplazamiento');
xlabel('Tiempo (seg.)');
ylabel('distancia');
grid on;

subplot(3, 2, 2);
hold on
plot(t,delta_p);
plot(t,delta_p_so);
hold off
legend({'Con observador','Sin observador'})
title('Velocidad');
xlabel('Tiempo (seg.)');
ylabel('Velocidad (m/s)');
grid on;

subplot(3, 2, 3);
hold on
plot(t,theta*(180/pi));
plot(t,theta_so*(180/pi));
hold off
legend({'Con observador','Sin observador'})
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 2, 4);
hold on
plot(t,theta_p);
plot(t,theta_p_so);
hold off
legend({'Con observador','Sin observador'})
title('Velocidad angular \omega_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 1, 3);
hold on
plot(t,u);
plot(t,u_so);
hold off
legend({'Con observador','Sin observador'})
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;


