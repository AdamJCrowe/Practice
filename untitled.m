%%%This programme is for  Worked Example 2 of lecture 1
clear all;
close all;
clc;

L=1; lambda=400; Ta=300; Tb=320; %physical constants
Sc=5000; Sp=-100; %source
n=21; %number of nodes

x0=linspace(0,L,n); %create grid of nodes
dx=L/(n-1); %distance between nodes
Dx=dx; %Dirichlet boundary conditions; we do not need the width of the boundary CV

A=zeros(n,n); B=zeros(n,1); %create matricies

A(1,1)=1; B(1)=Ta; %first node
A(n,n)=1; B(n)=Tb; %last node

for i=2:n-1 %fill matricies using loop from second node up to second from last node
    A(i,i-1)=-lambda/dx;
    A(i,i)=2*lambda/dx-Sp*Dx;
    A(i,i+1)=-lambda/dx;
    B(i)=Sc*Dx;
end

T=A\B; %find temprature profile using backlash
figure('color','w'); plot(x0,T,'rs'); grid on; %plot distance along volume against temprature at that node
xlabel('x [m]'); ylabel('T [K]'); hold on



%%%Tests to find accuarcy:

%Compare with theoretical solution, first some constants must be defined
%for the anaylitcal solution:
mu1=sqrt(abs(Sp)/lambda); mu2=-mu1;
c1=(Tb-(Sc/Sp+Ta)*exp(mu2*L)+Sc/Sp)/(exp(mu1*L)-exp(mu2*L)); c2=Ta+Sc/Sp-c1;
Tteo=c1*exp(mu1*x0)+c2*exp(mu2*x0)-Sc/Sp; %this is analyitcal solution
plot(x0,Tteo,'k-');
legend('Backlash','Exact');
error=mean(abs(T-Tteo')) %calulating error, ' next to Tteo transposes matrix

%Residual error check
residual=sum(abs(B-A*T))/sum(abs(diag(A).*T))

%Energy balance check
num=0; den=0;
for i=2:n-1
    num=num+abs(-lambda*(T(i)-T(i-1))/dx+lambda*(T(i+1)-T(i))/dx+(Sc+Sp*T(i))*dx);
    den=den+abs(lambda*(T(i)-T(i-1))/dx);
end
unbalance=num/den*100