function aicol(n,nf,np)
dt=0.1;
%defining variables
% T_c=Temperature of Glass cover
% T_a=Temperature of air gap
% T_ab=Temperature of absorber
% T_f=Temperature of fluid
% T_it=Temperature of top surface of insulation
% T_ib=Temperature of bottom surface of insulation
% T_am=Ambient Temperature
Pr=0.7;
T_am=293;
T_c=ones(n,1)*293;
T_a=ones(n,1)*293; 
T_ab=ones(n,1)*293;     
T_fi1=ones(n,1)*293; 
T_fi2=ones(n,1)*293;    
T_fi3=ones(n,1)*293;
T_fi4=ones(n,1)*293;
T_fi5=ones(n,1)*293;
T_f=ones(n,1)*293;  
T_it=ones(n,1)*293;
T_ib=ones(n,1)*293;
%Defining Temporary variable
hc_am=zeros(n,1); h_r1=zeros(n,1); hc_a=zeros(n,1); ha_ab=zeros(n,1); hab_f=zeros(n,1); hfi_f=zeros(n,1); hf_it=zeros(n,1); hi_am=zeros(n,1);
hc1=zeros(n,1); hc2=zeros(n,1); hc3=zeros(n,1); hc4=zeros(n,1); RR1=zeros(n,1); RR2=zeros(n,1); 
Nu_a=zeros(n,1);Nu_a1=zeros(n,1); Ra_1=zeros(n,1); RR3=zeros(n,1); RR4=zeros(n,1); Ra_2=zeros(n,1);
B=zeros(n,1); C=zeros(n,1);D=zeros(n,1); E=zeros(n,1);F=zeros(n,1);G=zeros(n,1);H=zeros(n,1);I=zeros(n,1);J=zeros(n,1);K=zeros(n,1);L=zeros(n,1);M=zeros(n,1);
N=zeros(n,1);O1=zeros(n,1);O2=zeros(n,1);O3=zeros(n,1);O4=zeros(n,1);P1=zeros(n,1);P2=zeros(n,1);P3=zeros(n,1);P4=zeros(n,1);P5=zeros(n,1);
P6=zeros(n,1);Q1=zeros(n,1);Q2=zeros(n,1);R1=zeros(n,1);R2=zeros(n,1);S=zeros(n,1);T=zeros(n,1);U=zeros(n,1);V=zeros(n,1);W1=zeros(n,1);W2=zeros(n,1);
W3=zeros(n,1);X=zeros(n,1);Y=zeros(n,1);Z=zeros(n,1);
%defining thickness of various parts
d_c=0.003;              %thickness of glass cover
d_a=0.02;               %thickness of air gap
d_ab=0.002;            %thickness of absorber plate           
d_i=0.04;               %thickness of insulation
%defining density values
rho_c=2400;
rho_a=1.2;
rho_ab=2707;
rho_fi=2707;
rho_f=1.2;
rho_i=35;
%defining specific heat values
cp_c=800;
cp_a=1005.7;
cp_ab=896;
cp_fi=896;
%cp_f=1005.7;
cp_i=1030;
%air properties
nw_am=15.5*10^(-6);
mw_am=1.8462*10^(-5);
k_am=0.028;
%other properties of different material used
r=0.05;          %reflectivity of glass cover
tau=0.9;         %transmissivity of glass cover
r1=0.03;         %reflectivity of absorber
alpha1=0.97;
k_ab=200;            %conductivity of absorber
k_fi=200;           %conductivity of fin material
k_i=0.037;       %conductivity of insulation
theta=((pi*32)/180);
g=9.8;
e_c=0.85;
e_ab=0.97;
e_i=0.6;
sigma=5.670*10^(-8);
%collector specification
l=1.4;
b=1.0;
%defining fin variables
t=0.002;     %fin thickness
h=0.07;      %fin height
c=0.12-h;       %fin clearance
p=b/nf;           %fin pitch
dz=l/n;
delta=(4*l*b)/(2*(l+b));
%mass_dot=(4000/3600);  %(in kg/(m^(-2)*s^(-1)))
%v=2.0;
a=(p*(h+c)-t*h);
A=a*nf;
M_dot=0.03; %(in kg/s for total collector)
mass_dot=M_dot/A;   %(in kg/m2s for total collector)
%v_dot=v*a;
wp=2*p+2*c+4*h;
dh=(4*a)/wp;
m_dot=M_dot/nf;     %(in kg/s for one section)
v_dot=m_dot/rho_f;   %(in m3/s for one section)
v=v_dot/a;            %(velocity for one section)
%G_i=628.27;               %(solar insolation)
av=3.0;               %(wind velocity)
rey=(mass_dot*dh)/mw_am;
ff=0.059*(rey^(-0.2));
delta_P=(2*ff*((mass_dot)^2)*l)/dh;
%convergence=0;
T_c_old2=T_c; T_a_old2=T_a; T_ab_old2=T_ab;T_fi1_old2=T_fi1;T_fi2_old2=T_fi2; T_fi3_old2=T_fi3; T_fi4_old2=T_fi4; T_fi5_old2=T_fi5;
T_f_old2=T_f; T_it_old2=T_it; T_ib_old2=T_ib;
t_o_p=(l*a)/v_dot;
tot=round(t_o_p/dt);
T_tot=np*tot;
fopen('temp2.txt');
fileID = fopen('temp2.txt','wt');
fprintf(fileID,'number of nodes = %6.1f \n flowrate = %6.3f \n ', n, M_dot);
fprintf(fileID,'T_c  T_a  T_ab   T_fi1   T_fi2   T_fi3   T_fi4   T_fi5   T_f   T_it   T_ib\n');
%calculating Temperatures
for k=1:T_tot
    if k<=18000
        G_i=546.35;
    elseif k>18000 && k<36000
        G_i=628.27;
    elseif k>36000  && k<54000
        G_i=710.22;
    elseif k>54000 && k<72000
        G_i=764.85;
    elseif k>72000 && k<90000 
        G_i=819.484;
    else
        G_i=835.87;
    end
    convergence=0;
    T_c_old1=T_c_old2; T_a_old1=T_a_old2; T_ab_old1=T_ab_old2; T_fi1_old1=T_fi1_old2;T_fi2_old1=T_fi2_old2; T_fi3_old1=T_fi3_old2;
    T_fi4_old1=T_fi4_old2; T_fi5_old1=T_fi5_old2; T_f_old1=T_f_old2; T_it_old1=T_it_old2; T_ib_old1=T_ib_old2;
    
    while convergence<=11*n
        T_c_old2=T_c; T_a_old2=T_a; T_ab_old2=T_ab; T_fi1_old2=T_fi1;T_fi2_old2=T_fi2; T_fi3_old2=T_fi3; T_fi4_old2=T_fi4; T_fi5_old2=T_fi5;
        T_f_old2=T_f; T_it_old2=T_it; T_ib_old2=T_ib;  
        for j=1:n
           T_a1(j)=T_a(j)-273;
           mw_a(j)=(1.983+0.00184*(T_a1(j)-27))*10^(-5);
           k_a(j)=0.02624+0.0000758*(T_a1(j)-27);
           cp_a1(j)=1.0057+0.000066*(T_a1(j)-27);
           nw_a(j)=mw_a(j)/rho_a;
           cp_a(j)=cp_a1(j)*1000;
        end
        for j=1:n
        T_f1(j)=T_f(j)-273;
        mw_f(j)=(1.983+0.00184*(T_f1(j)-27))*10^(-5);
        k_f(j)=0.02624+0.0000758*(T_f1(j)-27);
        cp_f1(j)=1.0057+0.000066*(T_f1(j)-27);
        nw_f(j)=mw_f(j)/rho_f;
        cp_f(j)=cp_f1(j)*1000;
        end
for j=1:n
%calculation of hc_am
Re_am(j)=(delta*av)/nw_am;
Nu_c1(j)=0.86*(Re_am(j)^0.5)*(Pr^(1/3));
hc1(j)=(Nu_c1(j)*k_am)/delta;
T_sky=0.0552*(T_am^1.5);
hc2(j)=(sigma*e_c*((T_c(j))^4-T_sky^4))/(T_c(j)-T_am);
if T_c(j)<293.15
    hc_am(j)=hc1(j);
else
hc_am(j)=hc1(j)+hc2(j);
end
%calculation of h_r1
h_r1(j)=(sigma*((T_ab(j))^2+(T_c(j))^2)*(T_ab(j)+T_c(j)))/((1/e_ab)+(1/e_c)-1);
%calculation of hc_a
beta=(1/T_a(j));
Ra_1(j)=((g*beta*abs(T_ab(j)-T_c(j))*(d_a^3))/(nw_a(j)^2))*0.7;
RR1(j)=(1-(1708/(Ra_1(j)*cos(theta))));
RR2(j)=(((Ra_1(j)*cos(theta))/5830)^(1/3)-1);
if RR1(j)<0 && RR2(j) < 0
    Nu_a(j)=1;
elseif RR1(j)>0 && RR2(j)<0
    Nu_a(j)=1+1.44*(1-((1708*((sin(1.8*theta))^1.6))/(Ra_1(j)*cos(theta))))*RR1(j);
elseif RR2(j)>0 && RR1(j)<0
    Nu_a(j)=1+RR2(j);
else
    Nu_a(j)=1+1.44*(1-((1708*(sin(1.8*theta))^1.6)/(Ra_1(j)*cos(theta))))*RR1(j)+RR2(j);
end
hc_a(j)=(Nu_a(j)*k_a(j))/d_a;
%calculation of ha_ab
Ra_2(j)=((g*beta*abs(T_ab(j)-T_c(j))*(d_a^3))/(nw_a(j)^2))*0.7;
RR3(j)=(1-(1708/(Ra_2(j)*cos(theta))));
RR4(j)=(((Ra_2(j)*cos(theta))/5830)^(1/3)-1);
if RR3(j)<0 && RR4(j)<0
    Nu_a1(j)=1;
elseif RR3(j)>0 && RR4(j)<0
    Nu_a1(j)=1+1.44*(1-((1708*(sin(1.8*theta))^1.6)/(Ra_2(j)*cos(theta))))*RR3(j);
elseif RR4(j)>0 && RR3(j)<0
    Nu_a1(j)=1+RR2(j);
else
    Nu_a1(j)=1+1.44*(1-((1708*(sin(1.8*theta))^1.6)/(Ra_2(j)*cos(theta))))*RR3(j)+RR4(j);
end
ha_ab(j)=(Nu_a1(j)*k_a(j))/d_a;
%calculation of hab_f
Re_f1(j)=(dh*v)/nw_f(j);
lh=0.05*Re_f1(j)*dh*Pr;
if (j*dz)<=lh
    Nu_f1(j)=7.54+((0.03*(dh/l)*Re_f1(j)*Pr)/(1+(0.016*((dh/l)*Re_f1(j)*Pr)^(2/3))));
    hab_f(j)=(Nu_f1(j)*k_f(j))/dh;
else
Nu_f(j)=5.4+((0.00190*((Re_f1(j)*Pr*(dh/l))^1.71))/(1+0.00563*((Re_f1(j)*Pr*(dh/l))^1.17)));
hab_f(j)=(Nu_f(j)*k_f(j))/dh;
end
%calculation of hfi_f
if (j*dz)<=lh
    Nu_f1(j)=7.54+((0.03*(dh/l)*Re_f1(j)*Pr)/(1+(0.016*((dh/l)*Re_f1(j)*Pr)^(2/3))));
    hfi_f(j)=(Nu_f1(j)*k_f(j))/dh;
else
hfi_f(j)=(Nu_f(j)*k_f(j))/dh;
end
%calculation of hf_it
if(j*dz)<=lh
    Nu_f1(j)=7.54+((0.03*(dh/l)*Re_f1(j)*Pr)/(1+(0.016*((dh/l)*Re_f1(j)*Pr)^(2/3))));
    hf_it(j)=Nu_f1(j)*k_f(j)/dh;
else
hf_it(j)=(Nu_f(j)*k_f(j))/dh;
end
%calculation of hi_am
Re_am1(j)=(delta*av)/nw_am;
Nu_c2(j)=0.86*(Re_am1(j)^0.5)*(Pr^(1/3));
hc3(j)=(Nu_c2(j)*k_am)/delta;
T_sky1=0.0552*(T_am^1.5);
hc4(j)=(sigma*e_i*((T_ib(j))^4-T_sky1^4))/(T_ib(j)-T_am);
if T_ib(j)<293.15
hi_am(j)=hc3(j);
else
    hi_am(j)=hc3(j)+hc4(j);
end
%calculation of coefficients
B(j)=((1-r-tau+r1*tau)/(d_c*cp_c*rho_c));
C(j)=(hc_am(j)/(d_c*cp_c*rho_c));
D(j)=(hc_a(j)/(d_c*cp_c*rho_c));
E(j)=(h_r1(j)/(d_c*cp_c*rho_c));
S(j)=((1/dt)+D(j)+C(j)+E(j));
F(j)=(hc_a(j)/(d_a*cp_a(j)*rho_a));
G(j)=(ha_ab(j)/(d_a*cp_a(j)*rho_a));
T(j)=((1/dt)+F(j)+G(j));
H(j)=(h_r1(j)/(d_ab*cp_ab*rho_ab));
I(j)=(ha_ab(j)/(d_ab*cp_ab*rho_ab));
J(j)=((tau-r1*tau)/(d_ab*cp_ab*rho_ab));
K(j)=(((2*k_ab)/d_ab)+((10*k_fi)/h))*(t/(p*d_ab*rho_ab*cp_ab));
L(j)=((hab_f(j)*(p-t))/(p*d_ab*rho_ab*cp_ab));
U(j)=((1/dt)+H(j)+I(j)+K(j)+L(j));
O1(j)=(((2*k_ab)/d_ab)+((10*k_fi)/h))*(5/(h*rho_fi*cp_fi));
O2(j)=((25*k_fi)/((h^2)*rho_fi*cp_fi));
O3(j)=((2*hfi_f(j))/(t*rho_fi*cp_fi));
O4(j)=((5*hfi_f(j))/(h*rho_fi*cp_fi));
W1(j)=((1/dt)+O1(j)+O2(j)+O3(j));
W2(j)=((1/dt)+2*O2(j)+O3(j));
W3(j)=((1/dt)+O2(j)+O3(j)+O4(j));
P1(j)=((2*hfi_f(j)*h)/(5*(p*c+(p-t)*h)*rho_f*cp_f(j)));
P3(j)=((hfi_f(j)*t)/((p*c+(p-t)*h)*rho_f*cp_f(j)));
P4(j)=((hab_f(j)*(p-t))/((p*c+(p-t)*h)*rho_f*cp_f(j)));
P5(j)=(m_dot/((p*c+(p-t)*h)*dz*rho_f));
P6(j)=(p*hf_it(j))/((p*c+(p-t)*h)*rho_f*cp_f(j));
X(j)=((1/dt)+5*P1(j)+P3(j)+P4(j)+P5(j)+P6(j));
Q1(j)=((2*hf_it(j))/(d_i*cp_i*rho_i));
Q2(j)=((4*k_i)/((d_i^2)*rho_i*cp_i));
Y(j)=((1/dt)+Q1(j)+Q2(j));
R1(j)=((2*hi_am(j))/(d_i*cp_i*rho_i));
R2(j)=((4*k_i)/((d_i^2)*rho_i*cp_i));
Z(j)=((1/dt)+R1(j)+R2(j));
end
        T_c(1)=(((T_c_old1(1)/dt)+B(1)*G_i+C(1)*T_am+D(1)*T_a(1)+E(1)*T_ab(1))/S(1));
        T_a(1)=(((T_a_old1(1)/dt)+F(1)*T_c(1)+G(1)*T_ab(1))/T(1));
        T_ab(1)=(((T_ab_old1(1)/dt)+H(1)*T_c(1)+I(1)*T_a(1)+J(1)*G_i+K(1)*T_fi1(1)+L(1)*T_f(1))/U(1));
        T_fi1(1)=(((T_fi1_old1(1)/dt)+O1(1)*T_ab(1)+O2(1)*T_fi2(1)+(O3(1))*T_f(1))/W1(1));
        T_fi2(1)=(((T_fi2_old1(1)/dt)+O2(1)*T_fi1(1)+O2(1)*T_fi3(1)+(O3(1))*T_f(1))/W2(1));
        T_fi3(1)=(((T_fi3_old1(1)/dt)+O2(1)*T_fi2(1)+O2(1)*T_fi4(1)+(O3(1))*T_f(1))/W2(1));
        T_fi4(1)=(((T_fi4_old1(1)/dt)+O2(1)*T_fi3(1)+O2(1)*T_fi5(1)+(O3(1))*T_f(1))/W2(1));
        T_fi5(1)=(((T_fi5_old1(1)/dt)+O2(1)*T_fi4(1)+(O3(1)+O4(1))*T_f(1))/W3(1));
        T_f(1)=(((T_f_old1(1)/dt)+(P1(1))*(T_fi1(1)+T_fi2(1)+T_fi3(1)+T_fi4(1)+T_fi5(1))+P3(1)*T_fi5(1)+P4(1)*T_ab(1)+P5(1)*T_am+P6(1)*T_it(1))/X(1));
        T_it(1)=(((T_it_old1(1)/dt)+Q1(1)*T_f(1)+Q2(1)*T_ib(1))/Y(1));
        T_ib(1)=(((T_ib_old1(1)/dt)+R1(1)*T_am+R2(1)*T_it(1))/Z(1));
        for j=2:n
        T_c(j)=(((T_c_old1(j)/dt)+B(j)*G_i+C(j)*T_am+D(j)*T_a(j)+E(j)*T_ab(j))/S(j));
        T_a(j)=(((T_a_old1(j)/dt)+F(j)*T_c(j)+G(j)*T_ab(j))/T(j));
        T_ab(j)=(((T_ab_old1(j)/dt)+H(j)*T_c(j)+I(j)*T_a(j)+J(j)*G_i+K(j)*T_fi1(j)+L(j)*T_f(j))/U(j));
        T_fi1(j)=(((T_fi1_old1(j)/dt)+O1(j)*T_ab(j)+O2(j)*T_fi2(j)+(O3(j))*T_f(j))/W1(j));
        T_fi2(j)=(((T_fi2_old1(j)/dt)+O2(j)*T_fi1(j)+O2(j)*T_fi3(j)+(O3(j))*T_f(j))/W2(j));
        T_fi3(j)=(((T_fi3_old1(j)/dt)+O2(j)*T_fi2(j)+O2(j)*T_fi4(j)+(O3(j))*T_f(j))/W2(j));
        T_fi4(j)=(((T_fi4_old1(j)/dt)+O2(j)*T_fi3(j)+O2(j)*T_fi5(j)+(O3(j))*T_f(j))/W2(j));
        T_fi5(j)=(((T_fi5_old1(j)/dt)+O2(j)*T_fi4(j)+(O3(j)+O4(j))*T_f(j))/W3(j));
        T_f(j)=(((T_f_old1(j)/dt)+(P1(j))*(T_fi1(j)+T_fi2(j)+T_fi3(j)+T_fi4(j)+T_fi5(j))+P3(j)*T_fi5(j)+P4(j)*T_ab(j)+P5(j)*T_f(j-1)+P6(j)*T_it(j))/X(j));
        T_it(j)=(((T_it_old1(j)/dt)+Q1(j)*T_f(j)+Q2(j)*T_ib(j))/Y(j));
        T_ib(j)=(((T_ib_old1(j)/dt)+R1(j)*T_am+R2(j)*T_it(j))/Z(j));
        end
        error=zeros(n,1);
        ccc=0;
        for j=1:n
            if ccc<=0
            error(1)=abs(T_c(j)-T_c_old2(j))/T_c(j);
            error(2)=abs(T_a(j)-T_a_old2(j))/T_a(j);
            error(3)=abs(T_ab(j)-T_ab_old2(j))/T_ab(j);
            error(4)=abs(T_fi1(j)-T_fi1_old2(j))/T_fi1(j);
            error(5)=abs(T_fi2(j)-T_fi2_old2(j))/T_fi2(j);
            error(6)=abs(T_fi3(j)-T_fi3_old2(j))/T_fi3(j);
            error(7)=abs(T_fi4(j)-T_fi4_old2(j))/T_fi4(j);
            error(8)=abs(T_fi5(j)-T_fi5_old2(j))/T_fi5(j);
            error(9)=abs(T_f(j)-T_f_old2(j))/T_f(j);
            error(10)=abs(T_it(j)-T_it_old2(j))/T_it(j);
            error(11)=abs(T_ib(j)-T_ib_old2(j))/T_ib(j);
            for i=1:11
            if error(i)<=10^(-8)
                convergence=convergence+1;
            else
                ccc=1;
            end
            end
            end
        end
    end
    fprintf('%f\n',k)
    for i=1:n
        fprintf('%f T_c=%6.1f T_a=%6.1f T_ab=%6.1f T_fi1=%6.1f T_fi2=%6.1f T_fi3=%6.1f T_fi4=%6.1f T_fi5=%6.1f T_f=%6.1f T_it=%6.1f T_ib=%6.1f\n',i ,T_c(i), T_a(i), T_ab(i), T_fi1(i), T_fi2(i), T_fi3(i), T_fi4(i), T_fi5(i), T_f(i), T_it(i), T_ib(i)) 
    end
    if rem(k,tot)==0
        fprintf(fileID,'%f\n',k);
       for i=1:n
        fprintf(fileID,'%f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',i ,T_c(i), T_a(i), T_ab(i), T_fi1(i), T_fi2(i), T_fi3(i), T_fi4(i), T_fi5(i), T_f(i), T_it(i), T_ib(i));
       end
       fprintf(fileID,'Insolation=%f ',G_i);
       eff=(cp_f(j)*M_dot*(T_f(n)-T_am))/(G_i*tau*alpha1*l*b)*100;
       fprintf(fileID,'effectiveness in percentage=%f\n',eff);
      fprintf(fileID,'-------------------------------------------\n'); 
    end
end
if rem(k,18000)==0
for i=1:n
fprintf(fileID,'heat tr coeff=%f \n',hfi_f(i));
end
end
fprintf(fileID,'Insolation=%f ',G_i);
fprintf(fileID,'h=%f ',h);
fprintf(fileID,'t=%f ',t);
fprintf(fileID,'p=%f ',p);
fprintf(fileID,'A=%f ',A);
fprintf(fileID,'M_dot=%f ',M_dot);
fprintf(fileID,'v=%f ',v);
fprintf(fileID,'mass_dot=%f ',mass_dot);
fprintf(fileID,'hydraulic diameter=%f ',dh);
fprintf(fileID,'reynolds no=%f ',rey);
fprintf(fileID,'friction factor=%f ',ff);
fprintf(fileID,'Pr drop=%f ',delta_P);
%fprintf(fileID,'time of operation=%f\n',t_o_p1);
fclose(fileID);