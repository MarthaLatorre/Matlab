close all; clear all; clc;
%% datos
v=3.5109e-6; PE=9.81; %% PE=9,81KN/m3
L = 208; L1=0.2*L; L2=0.25*L1; L3= 0.35*L;
E=4.6e-5
Kc_s=0.3; Kc_in=0.7; K_T=0.2; Ksal = 1; Kv= 0.05

g=9.81;

%% TUBERIA B
DB = 0.2429

%% tubería b
Reb=7.2e5; Erb=2.3e-4; Lb= L1+L+L3

syms F
eqn=-(1/(sqrt(F)))-2*log10(((Erb)/3.7)+(2.51/(Reb.*(sqrt(F)))))==0;
fb=solve(eqn,F);
fb=double(vpa(fb));

Db=E/Erb; Vb=(Reb*v)/Db; Qb=(Vb*pi()*(Db^2)/4);
Kb= (2*K_T)+(2*Kc_in);
hb=((fb*Lb/Db)+ Kb)*((Vb^2)/(2*g))

%Perdidas en las tuberias paralelas
hl=hb;

%% tubería rama  a
Da=0.0972; La=L1+L+L2
Ka=(2*K_T)+(2*Kc_s);

fa=0.03

syms Va
eqn= ((fa*La/Da)+Ka)*((Va^2)/(2*g)) - hl == 0 ; 
Va=solve(eqn,Va, 'Real', true); r=subexpr(Va(2)); 
Va=abs(double(vpa(r))); 
error=1;i=0;
while error>=0.01
 i=i+1;
M(i)=i;
% f(i)=f;

syms Va
eqn= ((fa*La/Da)+Ka)*((Va^2)/(2*g)) - hl == 0 ; 
Va1=solve(eqn(i),Va); 
r(i)=subexpr(Va1(2)); 
Va3(i)=(eval(r(i)));


Era(i)=E/Da;
Rea(i)=(Va3(i)*Da)/v;

if Rea(i)>= 4000
syms F
eqn(i)=-(1./(sqrt(F)))-2*log10(((Era(i))./3.7)+(2.51./(Rea(i).*(sqrt(F)))))==0;
fn=solve(eqn(i),F);
fna(i)=double(vpa(fn));

%f_new(i)=eval(fn(i));
error=abs((fna(i)-fa(i))/fna(i));
fa(i+1)=fna(i);



% % % FLUJO LAMINAR
else
   f_new(i)=64/Re(i)
   error=abs((f_new(i)-f(i))/f_new(i));
    f(i+1)=f_new(i);
 end 
end
fa 
Rea 
Va

 Qa=Va3(i)*pi()*(Da^2)*0.25
 

 T=[fa(i)'  Era(i)' Va3(i)'  Rea(i)' fna(i)]
fig = figure('Name','Iteraciones para hallar la velocidad en la rama supreior','Position',[400 100 500 500]);
cnames = {'f sup a','Rugosidad rel','Velocidad a','Reynolds a ','f nuevo a'};
t1 = uitable(fig,'Data',T,'ColumnName',cnames,'RowName',M','Position',[10 1 1000 500]);


%% TUBER�?A B
Q = Qa +Qb
VB = 4*Q/(pi()*DB^2);
ReB= VB*DB/v;
ErB=E/DB;

syms F
eqn=-(1/(sqrt(F)))-2*log10(((ErB)/3.7)+(2.5/(ReB*(sqrt(F)))))==0;
fn=solve(eqn,F);
fnB=double(vpa(fn))
hB=((fnB*L3/DB)+Kv)*(VB^2/(2*g))
hL12=hb+hB

PB= (152.4 - hL12 - (VB^2/(2*g)))*PE %%KPa

%% TUBERIA C
KC = 2*Kc_s + Ksal; LC =30+365.76; %% 365.76m son 1200ft
hC = ((fnB*LC/DB)+KC)*((VB^2)/(2*g))
h_A = 365.76 - (PB/PE) + ((VB^2)/(2*g)) + hC

%%PTENCIA AGREGADA
P_A= h_A*PE*Q %%KiloWatt




