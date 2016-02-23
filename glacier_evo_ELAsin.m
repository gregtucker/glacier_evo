%% glacier evolution finite difference diffusion solution
% Computer Modeling, Feb 2016, Week 5, JSB

clear all
figure(2)
clf

%% Initialize 

% Constants and variables
usl = 0.1; % sliding rate m/yr
A = 2.1e-16;  % this is available in any of the papers i sent along
Rho_i = 917; % ice  density
g = 9.81; % gravity

% arrays and variables

% time array
dt = .001; % time step 
tmax = 2 * 100; % max time , years
t = 0:dt:tmax; % time array

imax = length(t);

% x distance array
dx = 100; % horizontal distance step (m)
xmax = 10*1000; % horizontal distance max (m)
x = 0:dx:xmax; % horizontal distance array

% Setting up initial H
Hnaught = 0; % initial height of glacier
H = Hnaught * ones(size(x)); % initial H(x)

% Initial topography 
zbmax = 3000; % max bedrock
S = 0.15; % Slope of initial bedrock
zb = zbmax-(S*x); % initial bedrock topo
z = zb+H; % Initial valley topo

% Setting up mass balance 'b'
ELAnaught = 2500; % initial ELA
amp = 0; % amplitude m
P = 100; % period years
gamma = 0.01; % m per year/m altitude
zmax = 3000; % max elevation (m)
b = gamma*((zmax-(S*x))-ELAnaught); % initial mass balance

% plot animation
n=10; %number of plots
tplot = tmax/n;
time = 0;

% Analytical soln x-terminus
xterm = ((zmax-ELAnaught)*2)/S;
xtermline = xterm * ones(size(z));

%% RUN
for i = 1:imax
    % Calculate ELA at each time step
ELA = ELAnaught + amp*sin((2*pi*t(i))/P);
    % Calculate slope at each point
dzdx = diff(z)/dx; % slope
    % mass balance
b = gamma*((zmax-(S*x))-ELA); % mass balance
    % H edge
dHdx = diff(H)/dx; % H at the edge...this is the difference in H, not H at the edge
Hedge =H(1:end-1)+(0.5*dHdx);
    % Q
Q = (usl*Hedge)+((A*((Rho_i*g*abs(dzdx)).^3)).*((Hedge.^5)/5)); %Q
% new defn of Hedge
Q = [0 Q 0]; % pad Q array
    % dQdx
dQdx = diff(Q)/dx; % should be right size now
    % Ice thickness
dHdt = b-dQdx; % change in ice thickness over time
H = H + (dHdt*dt);
H = max(H,0); % cannot have negative ice

% recalculate topo
z = zb+H; % glacier topography with fixed bedrock and changing ice height

%% PLOT
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(2)
plot(x/1000,zb,'r','linewidth',2) % plot surface over time
hold on
plot(x/1000,z,'linewidth',1) % plot bedrock over time 
hold on

   xlabel('Distance (km)','fontname','arial','fontsize',24) % x label
   ylabel('Elevation (m)','fontname','arial','fontsize',24) % y label
   set(gca,'fontsize',18,'fontname','arial') % axes number labels
   title(['Glacier evolution after ',num2str(t(i)),' years']) % title - accumulates model time
   axis([0 xmax/1000 1500 3200]) % hold axes constant
   pause(0.05)

end
end
plot(xtermline/1000,z,'b--', 'linewidth',2) % analytical x terminus - out of loop
    


