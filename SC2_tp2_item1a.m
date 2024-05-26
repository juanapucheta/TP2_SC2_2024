% JUANA PUCHETA NOGUERA Practico 2 - Sistemas de Control 2

% Sistema en variables de estado que controle el ángulo del motor, para
% consignas de pi/2 y –pi/2 cambiando cada 5 segundos
%(el TL es el descrito en la planilla de datos)
% comparando el desempeño con el obtenido
% con el PID digital del TP No1. Hallar el valor de integración Euler
% adecuado.
% Objetivo: acelerar la dinámica del controlador.
clc; clear all; close all; 
tabla=xlsread('Curvas_Medidas_Motor_2024.xlsx');
t_D=tabla(:,1); %Tiempo
y_D=tabla(:,2); %Velocidad angular
i=tabla(:,3); %Corriente de armadura
u=tabla(:,4); %Tension
Tl=tabla(:,5); %Torque
Laa= 6.2811e-04; J=2.2328e-09;Ra=28.131;Bm=0;Ki= 0.011619;Km=0.060530;

% Defino el modelo del motor como variables de estado
A= [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B= [1/Laa 0; 0 -1/J; 0 0];
C= [0 0 1];
D= [0];

% Para definir el timepo de integracion de Euler busco el polo mas
% grande de la matriz A.
eig(A) % =  -2.2393e+04

% Calculo entonces el valor de integracion de Euler ADECUADO.
% te= (log(0.95)/-2.2393e+04)/10
% = 2.2906e-07 Redondeo a 3e-7

t_etapa=30e-7;

Q=diag([100 1e-3 1]);R=1e-1;
%Construcción del Hamiltoniano para el cálculo del controlador
Aa= A; Ba=B(:,1);
Ha=[Aa -Ba*inv(R)*Ba'; -Q -Aa'];
[n,va]=size(Ha); %6x6
[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:n/2,:); MX2=MX1X2(n/2+1:end,:);
% P=abs(MX2*inv(MX1));
P=real(MX2*inv(MX1));
Ka=inv(R)*Ba'*P;
%Fin cálculo del controlador-----------------------------------------------
disp('Controlador ampliado en ')
eig(Aa-Ba*Ka)

%Cálculo del Observador------------------------------------------
%Los valores del controlador de obtienen del K ampliado
%regulador a cero
T=.1;
e=0;ys=[];ypp(1)=0;yp(1)=0;y(1)=0;t=[];
Max_T=T/t_etapa;
y(1)=0;
X=[0;0;10]; %[corriente; velocidad angular; tita]
ref=1;psi=0;
u=0;x_hat=[0;0];y_o(1)=0;acc(1)=0;
for ii=1:Max_T
    u=-Ka*[X];
    y1=C*X; 
    psi_p=ref-y1;
    psi=psi+psi_p*t_etapa;
    X=modmotorpruebatita(t_etapa, X, [u, 0]);

    t(ii)=ii*t_etapa;
    ys(ii)= X(3);
    acc(ii)=u;
end
figure(1); 
subplot(3,1,1);plot(t,ys,'k');title('Salida')
subplot(3,1,3);plot(t,acc,'k');title('Acción de control');xlabel('Tiempo [seg.]')
hold on; 

%Torque del TP1 con Ts normalizado 
th=0:t_etapa:t_D(end);
ent_tm=zeros(size(th));
for ii=1:length(th)

  if th(ii)>0.1863

  ent_tm(ii)=1e-3;
  end

  if th(ii)>0.3372

  ent_tm(ii)=0;
  end
  if th(ii)>0.4866

  ent_tm(ii)=1e-3;
  end
end

%Referencia distinta de cero. 
Aa=[A,[0;0; 0];-C,0]; %Matriz A ampliada
Ba=[Ba;0]; %Matriz B ampliada 
Mat_Aa=Aa;Mat_Ba=Ba;

%Cálculo del LQR
Q=diag([1 1e-3 1 1e3]);R=1.5e-2;
%Construcción del Hamiltoniano para el cálculo del controlador
%Aa= A; Ba=B(:,1);
Ha=[Aa -Ba*inv(R)*Ba'; -Q -Aa'];
[n,va]=size(Ha); %6x6
[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:n/2,:); MX2=MX1X2(n/2+1:end,:);
% P=abs(MX2*inv(MX1));
P=real(MX2*inv(MX1));
Ka=inv(R)*Ba'*P;
%Fin cálculo del controlador-----------------------------------------------
disp('Controlador ampliado en ')
eig(Aa-Ba*Ka)


%Cálculo del Observador con referencia pi/2------------------------------------------
%Los valores del controlador de obtienen del K ampliado
%regulador a cero
T=.6;
e=0;ys=[];yi=[];ypp(1)=0;yp(1)=0;y(1)=0;t=[];
Max_T=T/t_etapa;
y(1)=0;
X=[0;0;0]; %[corriente; velocidad angular; tita]
ref=pi/2;psi=0;
u=0;x_hat=[0;0];y_o(1)=0;acc(1)=0;
for ii=1:Max_T
    u=-Ka*[X; psi];
    y1=C*X; 
    psi_p=ref-y1;
    psi=psi+psi_p*t_etapa;
    X=modmotorpruebatita(t_etapa, X, [u, ent_tm(ii)]);

    t(ii)=ii*t_etapa;
    yi(ii)= X(1);
    ys(ii)= X(3); 
    acc(ii)=u;
end
figure(2);
subplot(3,1,1);plot(t,ys,'k');title('Salida')
subplot(3,1,2);plot(t,yi,'k');title('Corriente')
subplot(3,1,3);plot(t,acc,'k');title('Acción de control');xlabel('Tiempo [seg.]')
hold on;

%Cálculo del Observador con ref -pi/2
%Los valores del controlador de obtienen del K ampliado
%regulador a cero
T=.6;
e=0;ys=[];yi=[];ypp(1)=0;yp(1)=0;y(1)=0;t=[];
Max_T=T/t_etapa;
y(1)=0;
%X=[0;0;0]; %[corriente; velocidad angular; tita]
ref=-pi/2;
u=0;x_hat=[0;0];y_o(1)=0;acc(1)=0;
for ii=1:Max_T
    u=-Ka*[X; psi];
    y1=C*X; 
    psi_p=ref-y1;
    psi=psi+psi_p*t_etapa;
    X=modmotorpruebatita(t_etapa, X, [u, ent_tm(ii)]);

    t(ii)=ii*t_etapa;
    yi(ii)= X(1);
    ys(ii)= X(3); 
    acc(ii)=u;
end
figure(3);
subplot(3,1,1);plot(t,ys,'k');title('Salida')
subplot(3,1,2);plot(t,yi,'k');title('Corriente')
subplot(3,1,3);plot(t,acc,'k');title('Acción de control');xlabel('Tiempo [seg.]')
hold on; 



