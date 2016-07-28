 clear
 close all
 clc
 global ms1 ms2 mp1 mp2 Is1 Is2 R1 R2 lpend1 lpend2 G zeta1 zeta2
%Top hoop parameters
ms1 = 1.2;          %mass of the spherical shell
%mass of the internal mechanism
mp1 = 1.0;       %mass of the pendulum
Is1 = .18;         %moment of inertia of the spherical shell
R1 = .15;           %radius of the outer spherical shell
lpend1 = .12;       %lenght of the pendulum
G =-9.8;
%Bottom hoop parameters
ms2 = 4*1.2;          %mass of the spherical shell
mp2 = 2.0;       %mass of the pendulum
Is2 = 2*.18;         %moment of inertia of the spherical shell
R2 = 2*.15;           %radius of the outer spherical shell
lpend2 = 2*.12;       %lenght of the pendulum
zeta1 =10; % friction coefficients
zeta2 =10;
global q;
%Initial position conditions

theta1 = 0*pi/180;
theta2 = 0*pi/180;

%Ode45 settings (for better results)
options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4 ,1e-4,1e-4,1e-4,1e-4]); %bottom hoop
%X = [theta1;theta2;theta1_dot;theta2_dot;phi1;phi2;phi1_dot;phi2_dot;]
x0 = [theta1-theta2*R2/R1; theta2; 0; 0; -20*pi/180; 20*pi/180; 0; 0]; % initial angles
%Simulate the system
tf = 10; % final time
tstart =0;% start time

[T,X] = ode45(@(t,X)ode2Hoops(t,X),linspace(tstart,tf,(tf-tstart)*100),x0);

plot(T,X)
legend 'theta1' 'theta2' 'theta1_dot' 'theta2_dot' 'phi1' 'phi2' 'phi1_dot' 'phi2_dot';


simtimescale = 4; % by what factor to scale down the simulation playback wrt ODE solution

animeFigure = figure();
animeFigure.Name = 'Simulation';
i = 1;
j = 1;
mov(1:round(length(T)/simtimescale)) = struct('cdata', [],'colormap', []);
while i <=length(T)    %Display every frame
    showHoopOnOtherHoop(animeFigure,[0,0],R1,X(i,1), X(i,5),lpend1,[0,0],R2,X(i,2), X(i,6),lpend2); % showSate plots the frame
    %pause(1/(tf*100))
    % capture the frame in the animation in video object
    mov(j)=getframe(animeFigure);
    j= j+1;
    i= i+simtimescale;
end

movie2avi(mov, 'Run2.avi', 'compression', 'none','quality',5)