close all; clear all; clc;
% anuncio=menu('Seleccione un elemento','Serie clase I','Serie clase II','Serie clase III','Paralelo');

Titulo = 'Dm, Tubería en serie o paralelo';
        NV = {'I n g r e s e   e l   C a u d a l','I n g r e s e   l a   r u g o s i d a d   d e l   m a t e r i a l ',...
            'I n g r e s e   l a   l o n g i t u d   d e   l a   t u b e r í a','I n g r e s e   l a   v i s c o s i d a d  c i n e m a t i c a                       '...
            'I n g r e s e   e l   p e s o   e s p e c í f i c o   d e   l a   s u s t a n c i a ','I n g r e s e   g r a v e d a d ',...
            'I n g r e s e  p o t e c i a  a g r e g a d a','I n g r e s e  p o t e c i a  r e t i r a d a'};
        NL = 1;
        %RD = {'0.8912','1.5e-04','75','9.15e-6','888888','32.2'};
        RD = {'0.07','4.6e-5','4.72','3.5137e-3','8.3607','9.81','0','0'};
        R = inputdlg(NV,Titulo,NL,RD);
        Q=str2double(R{1}); E=str2double(R{2}); L=str2double(R{3}); v=str2double(R{4}); PE=str2double(R{5});
        g=str2double(R{6}); hA=str2double(R{7}); hR=str2double(R{8});
        
        NV = {'Ingrese la P1','Ingrese la P2',...
            'Ingrese la z1','Ingrese la z2',...
            'Ingrese la V1','Ingrese la V2'};
        NL = 1;
        %RD = {'0','0','12','0','0','0'};
        RD = {'198.5','0','125','0','0','4.316'};
        R = inputdlg(NV,Titulo,NL,RD);
        P1=str2double(R{1}); P2=str2double(R{2}); z1=str2double(R{3}); z2=str2double(R{4});
        V1=str2double(R{5}); V2=str2double(R{6});
        
        NV = {'Ingrese un valor inicial para f'};  
        NL = 1; RD = {'0.06'}; R = inputdlg(NV,Titulo,NL,RD);
%         NL = 1; RD = {'0.06'}; R = inputdlg(NV,Titulo,NL,RD);
        f=str2double(R{1});
        
hl=(z1-z2)+((P1-P2)/PE)+((V1^2-V2^2)/(2*g))+hA-hR
syms D
% D=0.4011
%Pérdidas Mayores

hlM=(f*L*(8*Q^2))/(pi()^2*D^2*g);

%Pérdidas Menores
 NV = {'Ingrese la k1','Ingrese la k2',...
            'Ingrese la k3','Ingrese la k4',...
            'Ingrese la k5','Ingrese la k6'};
        NL = 1;
        %RD = {'0.5','10','0.3','0','0','0'};
        RD = {'1.2','0','0','0','0','0'};
        R = inputdlg(NV,Titulo,NL,RD);
        K1=str2double(R{1}); K2=str2double(R{2}); K3=str2double(R{3}); K4=str2double(R{4});
        K5=str2double(R{5}); K6=str2double(R{6});

Kt=K1+K2+K3+K4+K5+K6;

hlm=(Kt*V2^2)/(2*g);

% DESPEJANDO FACTOR DE FRICCIÓN

 %f=(z1*2*g*pi()^2*D^5)/(16*Q^2*L)- (Kt+1)*D/L
 Vc= ((8*(Q^2))/(pi()^2*(D^4)*g));

%  eqn=((2*z1*g*pi()^2*D^5)/(16*Q^2*L)- (Kt+1)*D/L - f == 0)
 
  eqn= ((L*f/D)+ Kt)*Vc - hl == 0 ; %% USE ((L*f(i)/D)+ Kt+1)*Vc(i) - hl == 0 SI CONOCE LA VELOCIDAD DE SALIDA
Dm=solve(eqn,D, 'MaxDegree', 5, 'PrincipalValue',true); %Principal Value toma la primera sol de D
Dm=double(vpa(Dm)); % Evalua D y convierte el sym en valor
error=1;i=0;
while error>=0.01
 i=i+1;
M(i)=i;
% f(i)=f;

syms D
Vc(i)= ((8*(Q^2))/(pi()^2*(D^4)*g));
eqn= ((L*f(i)/D)+ Kt)*Vc(i) - hl == 0; %% USE ((L*f(i)/D)+ Kt+1)*Vc(i) - hl == 0 SI CONOCE LA VELOCIDAD DE SALIDA
Dm(i)=solve(eqn,D, 'MaxDegree', 5, 'PrincipalValue',true); %Principal Value toma la primera sol de D
Dm(i)=(double(vpa(Dm(i))));


Er(i)=E/Dm(i);
V(i)= (4*Q)/(pi()*(Dm(i)^2));
Re(i)=(V(i)*Dm(i))/v;

if Re(i)>= 4000
syms F
eqn(i)=-(1./(sqrt(F)))-2*log10(((Er(i))./3.7)+(2.51./(Re(i).*(sqrt(F)))))==0;
fn=solve(eqn(i),F);
fn(i)=double(vpa(fn));
f_new(i)=eval(fn(i));
error=abs((fn(i)-f(i))/fn(i));
f(i+1)=fn(i);
%T=[f,D,Er,V,Re,fn];
% % % FLUJO LAMINAR
else
   f_new(i)=64/Re(i)
   error=abs((f_new(i)-f(i))/f_new(i));
    f(i+1)=f_new(i);
 end 
end

% f_new=eval(fn)
T=[f(1:i)' Dm' Er' V'  Re' f_new']
fig = figure('Name','Iteraciones para hallar el diametro','Position',[400 100 500 500]);
cnames = {'f sup','Diametro','Rugosidad rel','Velocidad 2','Reynolds','f nuevo'};
t1 = uitable(fig,'Data',T,'ColumnName',cnames,'RowName',M','Position',[10 1 1000 500]);
% M=['1'; '2'; '3'];
% f = figure('Name','Cantidad de Energía de Cada Pared (Mes-Hora)','Position',[100 100 1000 500]);
% cnames = {'f sup','Diametro','Rugosidad rel','Velocidad 2','Reynolds','f nuevo'};
% t1 = uitable(f,'Data',T,'ColumnName',cnames,'RowName',M,'Position',[10 1 1000 500]);