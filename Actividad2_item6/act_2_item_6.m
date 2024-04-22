clear; close all;

%--------------------------------------------
% Item 6 Trabajo Practico Control 2
%--------------------------------------------

X=[0; 0; 0];ii=0;t_etapa=1e-6;wRef=1;tF=0.6;

%Constantes del PID
Kp=10;Ki=1e3;Kd=50*1e-6;color_='r';
% Kp=1;Ki=0;Kd=0.0001;color_='k';
% Kp=10;Ki=0;Kd=0;color_='b';
% Kp=0;Ki=0;Kd=0;color_='b';

Ts=t_etapa;

t_delay = 0:Ts:0.0351;
t_1 = t_delay(end)+Ts:Ts:0.1869;
% t_1 = 0:Ts:0.1869;
t_2 = t_1(end)+Ts:Ts:0.3366;
t_3 = t_2(end)+Ts:Ts:0.4872;
t_4 = t_3(end)+Ts:Ts:0.6;

TL_v = [zeros(1,length([t_delay t_1])) 1e-3+1e-5*randn(1,length(t_2)) ...
          zeros(1,length(t_3)) 1e-3+1e-5*randn(1,length(t_4))];
      
% TL_v = [zeros(1,length(t_1)) 1e-3+1e-5*randn(1,length(t_2)) ...
%           zeros(1,length(t_3)) 1e-3+1e-5*randn(1,length(t_4))];
      
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

e=zeros(tF/t_etapa,1);u=0;

x1 = zeros(tF/t_etapa,1); 
x2 = zeros(tF/t_etapa,1);
x3 = zeros(tF/t_etapa,1);
acc = zeros(tF/t_etapa,1);

t_v = [t_delay t_1 t_2 t_3 t_4];
% t_v = [t_1 t_2 t_3 t_4];

for idx = 1:length(t_v)
    
 ii=ii+1;k=ii+2;
 
 TL = TL_v(idx); % perturbacion
 
 X=modmotor(t_etapa, X, u, TL);
 
 e(k)=wRef-X(3); %ERROR
 
 u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
 
 if idx < round(t_delay(end)/Ts)
     u = 0;
 end
 
 x1(ii)=X(1);%ia
 x2(ii)=X(2);%w
 x3(ii)=X(3);%theta
 acc(ii)=u;
end

% t=0:t_etapa:tF;

%% Plots

fz = 12;

figure

subplot(5,1,1);hold on;grid on;
plot(t_v,x1,color_);title('Salida , ia','Interpreter','latex','FontSize',fz);
ylabel('I[A]','Interpreter','latex','FontSize',fz);

subplot(5,1,2);hold on;grid on;
plot(t_v,x2,color_);title('Salida, Velocidad Angular','Interpreter','latex','FontSize',fz);
ylabel('w[rad/s]','Interpreter','latex','FontSize',fz);

subplot(5,1,3);hold on;grid on;
plot(t_v,x3,color_);title('Salida, Angulo','Interpreter','latex','FontSize',fz);
ylabel('theta[rad]','Interpreter','latex','FontSize',fz);

subplot(5,1,4);hold on;grid on;
plot(t_v,acc,color_);title('Entrada, Va','Interpreter','latex','FontSize',fz);
ylabel('V[V]','Interpreter','latex','FontSize',fz);

subplot(5,1,5);hold on;grid on;
plot(t_v, TL_v, color_);title('Entrada, TL','Interpreter','latex','FontSize',fz);
ylabel('TL[Nm]','Interpreter','latex','FontSize',fz);
xlabel('Tiempo [s]','Interpreter','latex','FontSize',fz);


% % Para verificar
% Laa=366e-6;
% J=5e-9;
% Ra=55.6;
% B=0;
% Ki=6.49e-3;
% Km=6.53e-3;
% num=[Ki]
% den=[Laa*J Ra*J+Laa*B Ra*B+Ki*Km ]; %wpp*Laa*J+wp*(Ra*J+Laa*B)+w*(Ra*B+Ki*Km)=Vq*Ki