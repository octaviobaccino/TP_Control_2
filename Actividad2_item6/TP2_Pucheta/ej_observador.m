clc;clear all;close all;


%Parametros de motor
Laa=5e-3;J=0.004;
Ra=0.2;Bm=0.005;
Ki=6.5e-5;
Km=0.055;

%Matrices de estado
A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1];%salida solo el angulo
D=[0];


%Matrices ampliadas
Aa=[A zeros(3,1); -C 0];
Ba=[B(:,1); 0];
Ca=[C 0];

%LQR
Q=1*diag([1 1/100 1/100 10000]);    R=100;
Q=1*diag([10 1/100 1/100 10000/2.5]);    R=10;

Kamp=lqr(Aa,Ba,Q,R);


%Observador--------------------------------------
Ao=A';
Bo=C';
Co=B';

Qo=1e0*diag([1 1/10 1/10]);    Ro=1;

Ko=lqr(Ao,Bo,Qo,Ro);
%------------------------------------------------

%Variables
tsim=50; 
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
    tLin(1,i)=1.15e-4;
% tLin(1,i)=0;
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
% figure(1)
% plot(t,tLin);
% title('Torque de entrada');

%condiciones iniciales
ia(1)=0;                %Corriente de armadura x1
theta(1)=0;             %Valocidad angular x2
omega(1)=0;             %Posicion angular x3

estados=[ia(1); omega(1); theta(1)];
estados_obs=[ia(1);omega(1);theta(1)];
    
Xop=[0 0 0]';
x=[ia(1) omega(1) theta(1)]';
xobs=[0 0 0]'; %inicializacion para el observador

psi(1)=0;
integracion(1)=psi(1);

for i=1:round(tsim/h)
    
    psi_p=ref(i)-Ca(1:3)*estados-Ca(4)*integracion;
    psi(i)=integracion+psi_p*h;
    
    u(i) = -Kamp(1:3)*estados_obs-Kamp(4)*psi(i);
    %Variables del sistema lineal
    ia(i)= x(1);
    omega(i)= x(2);
    theta(i)= x(3);
    
    x1_p=-Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x2_p=Ki*x(1)/J-Bm*x(2)/J-tLin(i)/J;
    x3_p=x(2);

    xp=[x1_p; x2_p; x3_p];
    x=x+h*xp;
    
    %------con Observador----------------------
  
    ia_0(i)= xobs(1);
    omega_0(i)= xobs(2);
    theta_0(i)= xobs(3);
    
    y_sal_o(i) = C * estados_obs;
    y_sal(i)   = Ca(1:3) * estados + Ca(4)*integracion;
    x_hat_p     = A*xobs+B*u(i)+Ko*(y_sal(i)-y_sal_o(i));
    xobs       = xobs + x_hat_p*h;
    %--------------------------------------------
    
    estados=[ia(i);omega(i);theta(i)];
    integracion=psi(i);
    estados_obs=[ia_0(i);omega_0(i);theta_0(i)];
    
end

%------Sistema sin observador---------------------------------------------

%condiciones iniciales
ia_so(1)=0;                %Corriente de armadura x1
theta_so(1)=0;             %Valocidad angular x2
omega_so(1)=0;             %Posicion angular x3

estados_so=[ia_so(1);omega_so(1); theta_so(1)];
x_so=[ia_so(1) omega_so(1) theta_so(1)]';

psi_so(1)=0;
integracion_so(1)=psi_so(1);

for i=1:round(tsim/h)
    
    psi_p_so=ref(i)-C*estados_so;
    psi_so(i)=integracion_so+psi_p_so*h;
    
    u_so(i) = -Kamp(1:3)*estados_so-Kamp(4)*psi_so(i);
    %Variables del sistema lineal
    ia_so(i)= x_so(1);
    omega_so(i)= x_so(2);
    theta_so(i)= x_so(3);
    
    x1_p=-Ra*x_so(1)/Laa-Km*x_so(2)/Laa+u_so(i)/Laa;
    x2_p=Ki*x_so(1)/J-Bm*x_so(2)/J-tLin(i)/J;
    x3_p=x_so(2);
    
    xp_so=[x1_p; x2_p; x3_p];
    x_so=x_so+h*xp_so;
    
    estados_so=[ia_so(i);omega_so(i);theta_so(i)];
    integracion_so=psi_so(i);
end

%-------------------------------------------------------------------------


subplot(2, 2, 1);
hold on
plot(t,ia,'b');
plot(t,ia_so,'r');
hold off
legend({'Con observador','Sin observador'},'Location','southeast')
title('Corriente de armadura i_a');
xlabel('Tiempo (seg.)');
ylabel('Corriente (A)');
grid on;

subplot(2, 2, 2);
hold on
plot(t,omega,'b');
plot(t,omega_so,'r');
hold off
legend({'Con observador','Sin observador'})
title('Velocidad angular \omega_r');
xlabel('Tiempo (seg.)');
ylabel('Velocidad angular (rad/s)');
grid on;

subplot(2, 2, 3);
hold on
plot(t,theta,'b');
plot(t,theta_so,'r');
plot(t,ref,'g--');
hold off
legend({'Con observador','Sin observador','Referencia'})
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(2, 2, 4);
hold on
plot(t,u,'b');
plot(t,u_so,'r');
hold off
legend({'Con observador','Sin observador'})
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;