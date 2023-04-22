clc;
clear all;
close all;

global N; 
global dx; 
global dt;
%%

% parameters
lambda = (3*1.6*0.28/1.28)*10^11;
mu = 3*1.6*(1-0.56)/(2*(1.28))*10^11;
kappa = 10^5;
rho = 7900;

% number of nodes
N = 101; %odd
n = N/2 - 1;
L = 1;

%%
% grid parameters
dx = L/(N-1);
dt = 1E-9;
T = 10^7;
nt = T/dt;
x = linspace(0,L,N);


%%
% initial solution matrix
UM1 = 2;
UM2 = 1;
% x1_temp = UM1*linspace(1,0,n+2);
% x2_temp = UM2*linspace(1,0,n+2);
% u_n_2 = [UM1*linspace(0,1,n+2), x1_temp(2:end)];
% u_n_1 = [UM2*linspace(0,1,n+2), x2_temp(2:end)];
u_n_1 = 0.01*sin(2*pi*x/L);
u_n_2 = 0.01*sin(2*pi*x/L);
%u_n = zeros(N,1);

% boundary conditions
u_n_1(1) = 0;
u_n_1(end) = 0;

u_n_2(1) = 0;
u_n_2(end) = 0;

%%
nt = 900000;
%% MAIN LOOP
   h1 = figure('units','normalized','outerposition',[0 0 1 1]);

   filename = "ST_model.gif";
   g_counter = 1;
for i = 1:nt
    Q = (mu+lambda/2)*d2UdX2(u_n_1).*(6*dUdX(u_n_1) + 2 + 3*((dUdX(u_n_1)).^2));
    Q = Q + kappa*(dUdX(u_n_1).*d3UdX2dT(u_n_2,u_n_1) + d3UdX2dT(u_n_2,u_n_1) + d2UdX2(u_n_1).*d2UdXdT(u_n_2,u_n_1));
    u_n = (dt^2)*Q/rho + 2*u_n_1 - u_n_2;
    
    I1 = 2.+(1.+dUdX(u_n)).^2;
    I2 = 1.+2*((1.+dUdX(u_n)).^2);
    I3 = (1.+dUdX(u_n)).^2;
    
    B11 = (1.+dUdX(u_n)).^2;
    % Using Constitutive relations
    dSigdI1 = (lambda/4)*((1.+dUdX(u_n)).^2-1) + (mu/2)*((1.+dUdX(u_n)).^2+1);
    dSigdI2 = -mu/2;
    dSigdI3 = 0;
    
    T11 = (2*((I2.*dSigdI2./sqrt(I3))+sqrt(I3).*dSigdI3))+((2*dSigdI1./sqrt(I3)).*B11)-((2*sqrt(I3).*dSigdI2)./B11);
    PV = T11.*(u_n-u_n_1)/(dt);
    
    
    % boundary conditions
    u_n(1) = 0;
    u_n(end) = 0;
    
    u_n_2 = u_n_1;
    u_n_1 = u_n;
    
    if mod(i,1000)==0
        subplot(2,1,1);
        plot(x,u_n);
        xlim([-0.05,1.05]);
        ylim([-0.02,0.02]);
        xlabel("x in m");
        ylabel("u in m");
        title("Saint-Venant Model");
        text(0.8,0.015,"t ="+ num2str(i*dt)+"sec");
        grid on;

        subplot(2,1,2); 
        plot(x,PV);
        xlim([-0.05,1.05]);
        ylim([-10^13,10^13]);
        grid on;
        xlabel("x in m");
        ylabel("Component of Poynting's Vector in Pa m/sec");
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if g_counter == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
          g_counter = 0;
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        pause(0.05);
        clf;
    end
end
%%










