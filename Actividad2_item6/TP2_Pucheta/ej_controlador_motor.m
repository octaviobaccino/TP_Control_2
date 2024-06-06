clc;clear all;close all;


%Parametros de motor
Laa=5e-3;J=0.004;
Ra=0.2;Bm=0.005;
Ki=6.5e-5;
Km=0.055;


%Matrices de estado
A = [-Ra/Laa -Km/Laa 0;Ki/J -Bm/J 0;0 1 0];
B = [1/Laa 0;0 -1/J;0 0]; %considerando el torque
C = [0 0 1];             %salida posicion

%Matrices ampliadas
Aa=[A zeros(3,1); -C 0];
Ba=[B(:,1); 0];
Ca=[C 0];


%LQR
Q=1*diag([1 1/100 1/100 10000]);    R=100;
Q=1*diag([10 1/100 1/100 10000/2.5]);    R=10;
Kamp=lqr(Aa,Ba,Q,R);

%Variables
tsim=100; 
h=1e-4; 
t=0:h:(tsim-h);

%Entradas del sistema
flag=1;
contador=0;
ref=zeros(1,round(tsim/h));
for i=(round(1/h)):1:(tsim/h)
    if flag==1 
    ref(1,i)=pi/2;
    contador=contador+1;
        if contador==round(20/h)
            flag=0;
            contador=0;
        end
    else 
    ref(1,i)=-pi/2;
    contador=contador+1;
             if contador==round(20/h)
                flag=1;
                contador=0;
            end
    end
   
end
% figure(1)
% plot(t,ref);
% title('Referencia');
flag=1;
contador=0;
tLin=zeros(1,round(tsim/h));
for i=(round(1/h)):1:(tsim/h)
    if flag==1 
    tLin(1,i)=1.15e-3;
% tLin(1,i)=10e-5;
    contador=contador+1;
        if contador==round(20/h)
            flag=0;
            contador=0;
        end
    else 
    tLin(1,i)=0;
    contador=contador+1;
             if contador==round(20/h)
                flag=1;
                contador=0;
            end
    end
   
end
figure(1)
plot(t,tLin);
title('Torque de entrada');


%condiciones iniciales
ia(1)=0;                %Corriente de armadura x1
theta(1)=0;             %Valocidad angular x2
omega(1)=0;             %Posicion angular x3

estados=[ia(1);
        omega(1);
        theta(1)];

    
Xop=[0 0 0]';
x=[ia(1) omega(1) theta(1)]';

psi(1)=0;
integracion(1)=psi(1);

for i=1:round(tsim/h)
    
    psi_p=ref(i)-C*estados;
    psi(i)=integracion+psi_p*h;
    
    u(i) = -Kamp(1:3)*estados-Kamp(4)*psi(i);
    %Variables del sistema lineal
    ia(i)= x(1);
    omega(i)= x(2);
    theta(i)= x(3);
    
    x1_p=-Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x2_p=Ki*x(1)/J-Bm*x(2)/J-tLin(i)/J;
    x3_p=x(2);
    
    xp=[x1_p; x2_p; x3_p];
    x=x+h*xp;
    
    estados=[ia(i);omega(i);theta(i)];
    integracion=psi(i);
end

figure(2)
subplot(3, 2, 1);
plot(t,ia);
title('Corriente de armadura i_a');
xlabel('Tiempo (seg.)');
ylabel('Corriente (A)');
grid on;

subplot(3, 2, 2);
plot(t,omega);
title('Velocidad angular \omega_r');
xlabel('Tiempo (seg.)');
ylabel('Velocidad angular (rad/s)');
grid on;

subplot(3, 2, 3);
hold on
plot(t,theta);
plot(t,ref);
hold off
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 2, 4);
hold on
plot(t,u);
hold off
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;

subplot(3, 2, 5);
hold on
plot(t,tLin);
hold off
title('Torque de entrada TL');
xlabel('Tiempo (seg.)');
ylabel('Nm');
grid on;



