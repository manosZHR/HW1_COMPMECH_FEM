clear; clc

D=0.1; k=1e-3; L=20;

nex=input('Give number of elements in x direction: ');

xfirst=0.;
xlast=20;
deltax_linear=(xlast-xfirst)/nex;
deltax_quad=deltax_linear/2;

x_linear=linspace(xfirst,xlast,nex+1);
x_quad=linspace(xfirst,xlast,2*nex+1);

%% Linear basis functions
u_a1=linear_A(nex); %Dirichlet boundary condition
u_b1=linear_B(nex); %Neumann boundary condition

u_a1_0=2*u_a1(1)-u_a1(2); u_b1_0=2*u_b1(1)-u_b1(2);
u_a1_end=2*u_a1(end)-u_a1(end-1); u_b1_end=2*u_b1(end)-u_b1(end-1);

%Augmented mesh with fake nodes
k=length(u_a1)+2;
u_a1_augm=zeros(k,1);
u_a1_augm(1)=u_a1_0;
u_a1_augm(k)=u_a1_end;
u_b1_augm(1)=u_b1_0;
u_b1_augm(k)=u_b1_end;
for i=2:k-1
    u_a1_augm(i) = u_a1(i-1);
    u_b1_augm(i) = u_b1(i-1) ;
end


for i=2:length(u_a1)+1
   
    dcdx_a1(i-1)= ( u_a1_augm(i+1)-u_a1_augm(i-1) )/(2*deltax_linear);
    dcdx_b1(i-1)= ( u_b1_augm(i+1)-u_b1_augm(i-1) )/(2*deltax_linear);
    
end
dcdx_b1(end)= 0; 

%% Quad basis functions
u_a2=quad_A(nex); %Dirichlet boundary condition
u_b2=quad_B(nex); %Neumann boundary condition

u_a2_0=2*u_a2(1)-u_a2(2); u_b2_0=2*u_b2(1)-u_b2(2);
u_a2_end=2*u_a2(end)-u_a2(end-1); u_b2_end=2*u_b2(end)-u_b2(end-1);

%Augmented mesh with fake nodes
k=length(u_a2)+2;
u_a2_augm=zeros(k,1);
u_a2_augm(1)=u_a2_0;
u_a2_augm(k)=u_a2_end;
u_b2_augm(1)=u_b2_0;
u_b2_augm(k)=u_b2_end;
for i=2:k-1
    u_a2_augm(i) = u_a2(i-1);
    u_b2_augm(i) = u_b2(i-1) ;
end


for i=2:length(u_a2)+1
   
    dcdx_a2(i-1)= ( u_a2_augm(i+1)-u_a2_augm(i-1) )/(2*deltax_quad);
    dcdx_b2(i-1)= ( u_b2_augm(i+1)-u_b2_augm(i-1) )/(2*deltax_quad);
    
end
dcdx_b2(end)= 0;

%% Plotting the results
figure(1)
plot(x_linear,u_a1','r')
hold on
plot(x_linear,u_b1','b')
title('Linear basis functions')
xlabel('x')
ylabel('C')
legend('dirichled','neumann')
hold off

figure(2)
plot(x_linear,dcdx_a1,'r')
hold on
plot(x_linear,dcdx_b1,'b')
title('Linear basis functions')
xlabel('x')
ylabel('dC/dx')
legend('dirichled','neumann')
hold off

figure(3)
plot(x_quad,u_a2','r')
hold on
plot(x_quad,u_b2','b')
title('Quad basis functions')
xlabel('x')
ylabel('C')
legend('dirichled','neumann')
hold off

figure(4)
plot(x_quad,dcdx_a2,'r')
hold on
plot(x_quad,dcdx_b2,'b')
title('Quad basis functions')
xlabel('x')
ylabel('dC/dx')
legend('dirichled','neumann')
hold off

figure (5)
plot(x_linear,u_a1','b')
hold on
plot(x_quad,u_a2','g')
title('Linear vs Quad basis functions')
xlabel('x')
ylabel('C')
legend('linear','quad')
hold off

figure (6)
plot(x_linear,dcdx_a1,'b')
hold on
plot(x_quad,dcdx_a2,'g')
title('Linear vs Quad basis functions')
xlabel('x')
ylabel('dCdx')
legend('linear','quad')
hold off