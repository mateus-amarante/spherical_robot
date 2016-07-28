
function []= RobustControllerFlatOneHoop(theta10,theta20,dtheta10,dtheta20,...
                              theta1f,theta2f,dtheta1f,dtheta2f,tf)
clc; close all;
% the nominal model parameter:
% parameters in the paper.
m1 = 1.2; 
m2 = 1.85; 
m3 = 2.05; 
R = .15; 
I1 = .018;
I2 = .0017;
I3 = .0006; 
g = 9.81;
l = 0.12;
eta = 0.03; 
% the nominal parameter vector b0 is
b1 = m1*R^2 + I1;
b2 = (m2 + m3)*R^2;
b3 =  m3*R*l;
b4 = m3*l^2 + I3 + I2;
b5 = eta;
b6 = m3*l;
b0 = [ b1; b2; b3; b4; b5; b6 ];
torque = [];%Output torque (it is filled in the sysODE function
%% Trajectory planning block
% Initial condition
x0 = [theta10, theta20, dtheta10, dtheta20, 0, 0];
%x0=[-0.5,0,-1,0.1];
x0e = [-0.7,0.5,-0.2,0,0]; % an error in the initial state.
xf=[theta1f, theta2f, dtheta1f, dtheta2f];
% The parameter for planned joint trajectory 1 and 2.
global a1 a2 % two polynomial trajectory for the robot joint
global tend
%tend = tf;
tend = 10;
nofigure=1;
a1 = planarArmTraj(theta10,dtheta10, theta1f, dtheta1f,tend, nofigure);
a2 = planarArmTraj(theta20,dtheta20, theta2f, dtheta2f,tend, nofigure);

options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]);
[T,X] = ode45(@(t,x)planarArmODEUncertain(t,x),[0 tf],x0,options);

statePlot = figure();
statePlot.Name ='State Variables';

hold on;
%plot state variables
plot(T,X(:,1:5));
legend('theta', 'phi', 'theta dot', 'phi dot', 'theta des');

figure('Name','theta1');
plot(T, X(:,1),'r-');
hold on
plot(T, X(:,5),'b-');
title('Theta_1 under Robust Control');
figure('Name','theta2');
plot(T, X(:,2),'r-');
hold on
plot(T, a2(1)+a2(2)*T+ a2(3)*T.^2+a2(4)*T.^3, 'b-');
title('Theta_2 under Robust Control');

%plot torque
%     sizetorque = size(torque)
%     TT = linspace(0,tf,sizetorque(2));
%     hold off;
%     torquePlot = figure();
%     torquePlot.Name = 'Torque Magnitude'
%     plot(TT,torque);
torquePlot = figure();

torquePlot.Name = 'Torque Magnitude';

%Hacking the ode to compute the torque values from its integral
UU = X(1:end,6)- circshift(X(1:end,6),[1,0]);
TT = T(1:end) - circshift(T(1:end),[1,0]);
UU_new = UU./TT;

plot(T(2:end),UU_new(2:end));
%axis([0 5 -5 5]);

    function [dx ] = planarArmODEUncertain(t,x)
        
        K = 0.001*eye(2);
        Lambda = 0.001*eye(2);
        vec_t = [1; t; t^2; t^3]; % cubic polynomials

        if( t >= tend )
            vec_tend = [1; tend; tend^2; tend^3]; % cubic polynomials
            theta_d = [a1'*vec_tend; a2'*vec_tend];
            dtheta_d = [0;0];
            ddtheta_d = [0;0];
        else
            theta_d = [a1'*vec_t; a2'*vec_t];
            %ref = [ref,theta_d];
            % compute the velocity and acceleration in both theta 1 and theta2.
            a1_vel = [a1(2), 2*a1(3), 3*a1(4), 0];
            a1_acc = [2*a1(3), 6*a1(4),0,0 ];
            a2_vel = [a2(2), 2*a2(3), 3*a2(4), 0];
            a2_acc = [2*a2(3), 6*a2(4),0,0 ];
            dtheta_d =[a1_vel*vec_t; a2_vel* vec_t];
            ddtheta_d =[a1_acc*vec_t; a2_acc* vec_t];
        end
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        
        %the true model
        m2t = m2+ 10*rand(1);% m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        I2t = I2 + (15/12)*rand(1);
        
        % the actual dynamic model of the system is characterized by M and
        % C
        m11 = (m1+m2+m3)*R^2 + I1;
        m12 = m3*R*l*cos(x(2));
        m21 = m3*R*l*cos(x(2));
        m22 = m3*l^2+I2+I3;
        Mmat = [ m11, m12; m21, m22];
        Cmat = [ eta, eta - m3*R*l*sin(x(2))*x(4);  
                 eta, eta];
        Gmat = [ 0; m3*g*l*sin(x(2))];
            
        invM = inv(Mmat);
        invMC = invM*Cmat;

        % compute the robust controller
        e = theta - theta_d;
        de = dtheta - dtheta_d;
        r = de + Lambda*e;
        v = dtheta_d - Lambda*e;
        a = ddtheta_d - Lambda*de;
        Y2= [ a(1), a(1), a(2)*cos(x(2)) - x(4)*v(2)*sin(x(2)), 0, v(1) + v(2), 0;
                 0,    0,          a(1)*cos(x(2)), a(2), v(1) + v(2), g*sin(x(2))];
        epsilon=4;
        rho =  10;%13.46; % see eq(27) in the paper.
        if norm(Y2'*r) > epsilon
            u = -rho* Y2'*r/norm(Y2'*r);
        else
            u= - rho* Y2'*r/epsilon;
        end
        tau = Y2*(b0 + u)- K*r;
        torque = [torque, tau(1,1)];
        dx = zeros(4,1);
        dx(1) = x(3);
        dx(2) = x(4);
        dx(3:4) = -invMC* x(3:4) - invM*Gmat + invM*[tau(1,1);tau(1,1)]; % because ddot theta = -M^{-1}(C * dot Theta) - M^{-1}*G + M^{-1}*tau
        dx(5) = dtheta_d(1,1);
        dx(6) = tau(1,1);
    end
%% Show animation and make a video of output

simtimescale = 2; % by what factor to scale down the simulation playback wrt ODE solution
animeFigure = figure();
animeFigure.Name = 'Simulation';

i = 1;
j = 1;

mov(1:round(length(T)/simtimescale)) = struct('cdata', [],'colormap', []);

while i <=length(T) %Display every frame
    showHoop(animeFigure,[0,0],R,X(i,1), X(i,2),l,0); % showSate plots the frame
    % capture the frame in the animation in video object

    mov(j)=getframe(animeFigure);
    j= j+1;
    i= i+simtimescale;

end
movie2avi(mov, 'Run1.avi', 'compression', 'none','quality',5)

end