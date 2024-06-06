clc;clear all;
m=.1;Fricc=0.1; long=0.6; g=9.8; M=.5;
%Condiciones iniciales
alfa(1)=.1; color='r';
alfa(1)=.5; color='g';
alfa(1)=.8; color='b';
ref=10;
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
%Versi�n linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0]
Mat_B=[0; 1/M; 0; -1/(long*M)]
Mat_C=[1 0 0 0]; %La salida monovariable es posici�n
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B];%Matriz Controlabilidad
%C�lculo del controlador por asignaci�n de polos
auto_val=eig(Mat_A);
c_ai=poly(auto_val);
Mat_W=[c_ai(4) c_ai(3) c_ai(2) 1;
 c_ai(3) c_ai(2) 1 0;
 c_ai(2) 1 0 0;
 1 0 0 0];
Mat_T=Mat_M*Mat_W;
A_controlable=inv(Mat_T)*Mat_A*Mat_T; %Verificaci�n de que T est� bien
%Ubicaci�n de los polos de lazo cerrado en mui:
mui(1)=-.7;mui(2)=-.7; mui(3)=-10 + 0.4i;mui(4)=conj(mui(3));
alfa_ia=poly(mui);
K=fliplr((alfa_ia(2:5)-c_ai(2:5)))*inv(Mat_T);
Gj=-inv(Mat_C*inv(Mat_A-Mat_B*K)*Mat_B);
eig(Mat_A-Mat_B*K) %Verifico que todos los polos est�n en el semiplano izquierdo
%lamdai=-20, por lo tanto el tiempo asociado es T1=50e-3, h puede estar
% entre T1/3 y T1/30, o sea entre 17e-3 y 1.7e-3
h=5e-3;
tiempo=(20/h);p_pp=0;tita_pp=0; t=0:h:tiempo*h;
omega=0:h:tiempo*h; alfa=0:h:tiempo*h; p=0:h:tiempo*h;
p_p=0:h:tiempo*h; u=linspace(0,0,tiempo+1);
while(i<(tiempo+1))
 estado=[p(i); p_p(i); alfa(i); omega(i)];
 u(i)=-K*estado+Gj*ref;
 p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
 tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
 p_p(i+1)=p_p(i)+h*p_pp;
 p(i+1)=p(i)+h*p_p(i);
 omega(i+1)=omega(i)+h*tita_pp;
 alfa(i+1)=alfa(i)+h*omega(i);
 y_sal(i)=Mat_C*estado;
 i=i+1;
end
figure(1);hold on;
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad �ngulo');hold on;
subplot(3,2,2);plot(t,alfa,color);grid on;title('�ngulo');hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posici�n carro');hold on;
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Acci�n de control');xlabel('Tiempo en Seg.');hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('�ngulo');ylabel('Velocidad angular');hold on;
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;