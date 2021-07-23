%% Dynamic simulation of marine robot

clear all;
close all;
clc;

%% Simulation parameters 

dt = 0.1; % step size
ts = 10 ; % simulation time
t = 0:dt:ts ; % time span

%% Initial Conditions

eta0 = [0;0;0]; % initial position and orientation of the vehicle
zeta0 = [0;0;0]; % initial vector of input commands.

eta(:,1) = eta0 ;
zeta(:,1) = zeta0 ;

%% Boat parameters

m = 10 ; % mass of vehicle is 10 kgms
Iz =  0.1 ; % Inertia of vehicle

xbc=0; ybc=0 ; % coordinates of mass center
rho = 1000;    % density of water  
cd = 0.4;      % coefficient of drag
dia = 0.2;     % diameter of hull in meter
area = pi*(dia^2)/4 ; % calculating area of hull
beta= 0.5;     % skin friction factor

%% State propagation

for i=1:length(t)
    u = zeta(1,i); v = zeta(2,i); r = zeta(3,i) ;
    
    %% Inertia matrix, N vector 
    
    D = [m,0,-m*ybc;
         0,m,m*xbc;
         -m*ybc,m*xbc,Iz+m*(xbc^2+ybc^2)];
    
    n_v = [-m*r*(v+xbc) + 2*0.5*rho*cd*area*(u^2) + beta*u;
            m*r*(u-ybc*r);
            m*r*(xbc*u-ybc*v)] ;
    
    %% input vector
    tau(:,i) = [1;0.5;0];
    
    %% Jacobian matrix
    
    psi = eta(3,i);
    J_eta = [cos(psi),-sin(psi),0;
             sin(psi),cos(psi), 0;
              0,       0,      1];
    
    zeta_dot(:,i) = inv(D)*(tau(:,i) - n_v) ;
    zeta(:,i+1) = zeta(:,i) + dt*zeta_dot(:,i);    % euler method of integration
    
    eta(:,i+1) = eta(:,i) + dt*( J_eta*(zeta(:,i) + dt*zeta_dot(:,i))) ;
      
end    

%% Animation 

l = 1.2 ; % length of boat
w = 1.2 ; % width of boat

bo_co = [-l/2,l/2,l/2,-l/2,-l/2;
         -w/2,-w/2,w/2,w/2,-w/2];   % boat coordinates

figure

for i=1:length(t)
    psi = eta(3,i);
    R_psi = [cos(psi),-sin(psi);
             sin(psi),cos(psi);];   % rotattion matrix
    
    v_pos = R_psi*bo_co ;
    fill(v_pos(1,:) + eta(1,i),v_pos(2,:) + eta(2,i),'g');
    hold on, grid on;
    axis([-1 3 -1 3]), axis square
    plot(eta(1,1:i),eta(2,1:i),'b-');
    legend('MR','Path'),set(gca,'fontsize',24)
    xlabel('x,[m]');
    ylabel('y,[m]');
    pause(0.1);
    hold off
         
end    