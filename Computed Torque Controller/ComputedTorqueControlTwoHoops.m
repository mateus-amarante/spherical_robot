function []= ComputedTorqueControlTwoHoops()

clear; clc; close all;

%Top hoop parameters
t_m1 = 1.2;          %mass of the spherical shell
t_m2 = 1.85;         %mass of the internal mechanism
t_m3 = 4*2.05;       %mass of the pendulum
t_I1 = .018;         %moment of inertia of the spherical shell
t_I2 = .0017;        %moment of inertia of internal mechanism
t_I3 = .0006;        %moment of inertia of the pendulum
t_R = .15;           %radius of the outer spherical shell
t_lpend = .12;       %lenght of the pendulum

%Bottom hoop parameters
b_m1 = 2*1.2;        %mass of the spherical shell
b_m2 = 2*1.85;       %mass of the internal mechanism
b_m3 = 2*4*2.05;     %mass of the pendulum
b_I1 = 2*.018;       %moment of inertia of the spherical shell
b_I2 = 2*.0017;      %moment of inertia of internal mechanism
b_I3 = 2*.0006;      %moment of inertia of the pendulum
b_R = t_R*2;         %radius of bottom hoop
b_lpend = 2*t_lpend; %lenght of the pendulum

gravity = 9.81;      %acceleration due to gravity
zeta = 0.03;         %viscous damping coeficient associated with the pendulum-sphere bearing
tf = 15;             %simulation end point

%Auxillary variables
%top hoop
t_b = t_m3*t_R*t_lpend;
t_d = t_m3*gravity*t_lpend;
t_Mt = t_m1+t_m2+t_m3;
%bottom hoop
b_b = b_m3*b_R*b_lpend;
b_d = b_m3*gravity*b_lpend;
b_Mt = b_m1+b_m2+b_m3;

%Diagonal members of M(1) //Notice they are constants
%top hoop
t_m11 = t_Mt*t_R^2 + t_I1;
t_m22 = t_m3*t_lpend^2 + t_I2 + t_I3;
%bottom hoop
b_m11 = b_Mt*b_R^2 + b_I1;
b_m22 = b_m3*b_lpend^2 + b_I2 + b_I3;

%Initial state of the state variables
%x0 = [t_theta; t_theta_dot; t_phi; t_phi_dot; t_theta_desired; t_torque;...
x0 = [0; 0; 0; 0; 0; 0;
     %b_theta; b_theta_dot; b_phi; b_phi_dot; b_theta_desired; b_torque]
      0; 0; 0; 0; 0; 0];

%Ode45 settings (for better results)
options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4 ,1e-4,1e-4,...%top hoop
                                         1e-4, 1e-4, 1e-4, 1e-4 ,1e-4,1e-4]); %bottom hoop

%Simulate the system 
[T,X] = ode45(@(t,x)sysODE(t,x),[0 tf],x0, options);

%Plot top hoop states separately
statePlot = figure();
statePlot.Name ='Top Sphere State Variables';
hold on;
%plot state variables
plot(T,X(:,1:end-7));
legend('theta', 'theta dot', 'phi', 'phi dot','desired q');
hold off;

%Plot bottom hoop states separately
statePlot = figure();
statePlot.Name ='Bottom sphere State Variables';
hold on;
%plot state variables
plot(T,X(:,7:end-1));
legend('theta', 'theta dot', 'phi', 'phi dot','desired q');
hold off;

%Plot torque
torquePlot = figure();
torquePlot.Name = 'Top sphere Torque Magnitude';
%hacking the ode to compute the torque values from its integral
UU = X(1:end,6)- circshift(X(1:end,6),[1,0]);
TT = T(1:end) - circshift(T(1:end),[1,0]);
UU_new = UU./TT;
plot(T(2:end),UU_new(2:end));

%Plot torque
torquePlot = figure();
torquePlot.Name = 'Bottom Torque Magnitude';
%hacking the ode to compute the torque values from its integral
UU = X(1:end,12)- circshift(X(1:end,12),[1,0]);
TT = T(1:end) - circshift(T(1:end),[1,0]);
UU_new = UU./TT;
plot(T(2:end),UU_new(2:end));

%%
    %Simulation funtion
    function [x_dot] = sysODE(t,x)
    %Describing the state variables of the top hoop
        t_theta = x(1);
        t_theta_dot = x(2);
        t_phi = x(3);
        t_phi_dot = x(4);
        
    %Describing the state variables of the bottom hoop
        b_theta = x(7);
        b_theta_dot = x(8);
        b_phi = x(9);
        b_phi_dot = x(10);
        
    %Calculating the slope of the tangent at the point of contact
    %of the top sphere on the bottom sphere
        gamma =   -t_theta * t_R/b_R;
    
    %Calculating the desired position, velocity and acceleration values for 
    %top and bottom hoop 
        
        %Computing a dynamic desired positon, velocity and acceleration for
        %the hoop of the top shell
        traj = getTrajecteoryPoint(t);
        %Computing a dynamic desired position for pendulum of the top shell
        t_phiDesired = (t_m1+t_m2+t_m3)*t_R*sin(gamma)/(t_m3*t_lpend);% term can be added for extra correction -  asin(traj(3)*I1/d) 
        
        %Setting the desired points for the top and bottom hoop
        %xd = [theta, theta_dot, theta_ddot, phi, phi_dot, phi_ddot]
        t_xd = [traj(1), traj(2), traj(3), t_phiDesired, 0, 0];
        b_xd = [0,0,0, 0,0,0];
        
    %Calculating the intertia, coriollis and gravity components of the 
    %top and bottom hoop 
         
        %Inertia Matrix - M(q) components 
        t_m12 = t_b*cos(t_phi-gamma);
        t_m21 = t_m12;
        b_m12 = b_b*cos(b_phi);
        b_m21 = b_m12;
        
        %Term used in N(q,q_dot)
        t_damp = zeta*(t_theta_dot+t_phi_dot);
        b_damp = zeta*(b_phi_dot+ b_theta_dot);
        %Coriolis and Gravity Matrix - N(q,q_dot)
        t_n1 = t_damp + t_Mt*t_R*gravity*sin(gamma) - t_b*sin(t_phi-gamma)*t_phi_dot^2;
        t_n2 = t_damp + t_d*sin(t_phi);
        b_n1 = b_damp - b_b*sin(b_phi)*b_phi_dot^2;
        b_n2 = b_damp + b_d*sin(b_phi);
        
    %Calculating the controller torque values to be applied on the 
    %the top and bottom hoop
     
        %Position and velocity gains of Shell
        t_kp1 = 15.5;
        t_kv1 = 15.5;
        b_kp1 = 15.5;
        b_kv1 = 15.5;
        
        %Position and velocity gains of Pendulum
        t_kp2 = 17;
        t_kv2 = 17;
        b_kp2 = 17;
        b_kv2 = 17;
        
        %Controller - Computed torque control
        t_e_ddot = [-t_kp1*(x(1)-t_xd(1)) -t_kp2*(x(3)-t_xd(4))] + [-t_kv1*(x(2)-t_xd(2)) -t_kv2*(x(4)-t_xd(5))];
        t_u1 = t_m11*t_e_ddot(1) + t_m12*t_e_ddot(2) + t_n1 + t_m11*t_xd(5) + t_m12*t_xd(6);
        b_e_ddot = [-b_kp1*(x(7)-b_xd(1)) -b_kp2*(x(9)-b_xd(4))] + [-b_kv1*(x(8)-b_xd(2)) -b_kv2*(x(10)-b_xd(5))];
        b_u1 = b_m11*b_e_ddot(1) + b_m12*b_e_ddot(2) + b_n1 + b_m11*b_xd(5) + b_m12*b_xd(6);
        %u2 = 21*e_ddot(1) + m22*e_ddot(2) + n2 + m21*xd(5) + m22*xd(6);
        
        t_u = t_u1;%select torque
        b_u = b_u1;%select torque
        
    %Calculating the f(x) and g(x) components for describing the state 
    %equations of the top and bottom hoop
        
        %Denominator of f(x) and g(x) ('g' is 'b' in the article)
        t_f_den = t_m12*t_m21 - t_m11*t_m22;
        t_g_den = -t_f_den;
        b_f_den = b_m12*b_m21 - b_m11*b_m22;
        b_g_den = -b_f_den;
        
        %Compute f(x) and g(x)
        t_f = [(t_m22*t_n1-t_m12*t_n2), (t_m11*t_n2-t_m21*t_n1)]./t_f_den;
        t_g = [t_m22-t_m12, t_m11-t_m21]./t_g_den;
        b_f = [(b_m22*b_n1-b_m12*b_n2), (b_m11*b_n2-b_m21*b_n1)]./b_f_den;
        b_g = [b_m22-b_m12, b_m11-b_m21]./b_g_den;
        
        %State equations which compute x_dot of the state variables x
        %for the top hoop
        x_dot(1) = x(2);                %theta_dot
        x_dot(2) = t_f(1) + t_g(1)*t_u; %theta_dot_dot
        x_dot(3) = x(4);                %phi_dot;
        x_dot(4) = t_f(2) + t_g(2)*t_u; %phi_dot_dot        
        x_dot(5) = t_xd(2);             %desired theta_dot
        x_dot(6) = t_u;                 %torque_dot
        
        %State equations which compute x_dot of the state variables x 
        %for the bottom hoop
        x_dot(7)  = x(8);                %theta_dot
        x_dot(8)  = b_f(1) + b_g(1)*b_u; %theta_dob_dot
        x_dot(9)  = x(10);               %phi_dot;
        x_dot(10) = b_f(2) + b_g(2)*b_u; %phi_dob_dot
        x_dot(11) = b_xd(2);             %desired theta_dot
        x_dot(12) = b_u;                 %torque_dot
        
        %Taking the transpose of the x_dot vector
        x_dot = x_dot';
        
    end
%%
% global variable to hold the latest value of computed theta position
q=0
    function [Xd] = getTrajecteoryPoint(t)
        tstart = 0;
        tend = tf*.75;
        %keep computing trajectory values for some time, stop few seconds
        %before the simulation ends
        if(t<tend)
            %quintic polynimial profile
            A = [1 tstart tstart^2 tstart^3 tstart^4 tstart^5
                0 1 2*tstart 3*tstart^2 4*tstart^3 5*tstart^4
                0 0 2 6*tstart 12*tstart^2 20*tstart^3
                1 tend tend^2 tend^3 tend^4 tend^5
                0 1 2*tend 3*tend^2 4*tend^3 5*tend^4
                0 0 2 6*tend 12*tend^2 20*tend^3];
            B = [0;0;0;25*pi/180;0;0];
            
            const = inv(A)*B;
            q = const(1) + const(2)*t + const(3)*t^2 + const(4)*t^3 +const(5)*t^4 + const(6)*t^5;
            q_dot = const(2) + 2*const(3)*t + 3*const(4)*t^2 + 4*const(5)*t^3 + 5*const(6)*t^4;
            q_ddot =  2*const(3) + 6*const(4)*t + 12*const(5)*t^2 + 20*const(6)*t^3;

            %Xd = [theta, theta_dot, theta_ddot, phi, phi_dot, phi_ddot]
            Xd = [q, q_dot, q_ddot, 0, 0, 0];
            
        else
            %just stabilize the system at where it was for rest of the few
            %second time.
            Xd = [q, 0, 0, 0, 0, 0];%theta, theta_dot, theta_ddot, phi, phi_dot, phi_ddot
        end
    end

%% 

%Show animation and make a video of output

videotimescale = 4; % by what factor to scale down the video wrt ODE solution
simtimescale = 40; % by what factor to scale down the simulation playback wrt ODE solution

animeFigure = figure();
animeFigure.Name = 'Simulation';
framesCompPerSec = length(T)/tf; % calculate frames per second
videoFrameRate = framesCompPerSec/videotimescale; % calculate frames rate of the video
animationframeTime = 1/framesCompPerSec; % calculate frames rate of the simulation playback

% create a video object
mov(1:(length(T)/videotimescale)) = struct('cdata', [],'colormap', []);
for i = 1:length(T)-1    %Display every frame
    showHoopOnOtherHoop(animeFigure,[0,0],t_R,X(i,1), X(i,3),t_lpend,[0,0],b_R,X(i,7), X(i,9),b_lpend); % showSate plots the frame
    pause(animationframeTime/simtimescale)
    
    if(mod(i,videotimescale) == 0)
        % capture the frame in the animation in video object
        %  mov(i/videotimescale)=getframe(animeFigure);
    end
end
% save the movie
% movie2avi(mov, 'Run2.avi', 'compression', 'none','fps',round(videoFrameRate*3),'quality',5);


end