function [X]=modmotor(t_etapa, xant, accion, perturbacion)
Laa=5.7586e-04; J=2e-9;Ra=28.13;B=0;Ki=0.0122;Km=0.0605;
Va=accion;
TL = perturbacion;
h=1e-6;
% omega= xant(1);
% wp= xant(2);
A_m = [-Ra/Laa -Km/Laa 0 ; Ki/J -B/J 0 ; 0 1 0];
B1_m = [1/Laa 0*-1/J 0]';
B2_m = [0*1/Laa -1/J 0]';
X = [xant(1) xant(2) xant(3)]';

ia = zeros(1,t_etapa/h);
w = zeros(1,t_etapa/h);
theta = zeros(1,t_etapa/h);

for ii=1:t_etapa/h
    
    ia(ii) = X(1);
    w(ii) = X(2);
    theta(ii) = X(3);
    
    Xp = A_m*X + B1_m.*Va + B2_m.*TL;
    X = X + Xp.*h;
%  wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
%  wp=wp+h*wpp;
%  omega = omega + h*wp;
end
% X=[omega,wp];
