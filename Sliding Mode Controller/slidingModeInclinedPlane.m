function []= slidingModeInclinedPlane()

    clear;
    clc;
    close all;

    m1 = 1.2; %mass of the spherical shell
    m2 = 1.85; %mass of the internal mechanism
    m3 = 2.05; %mass of the pendulum

    I1 = .018; %moment of inertia of the spherical shell
    I2 = .0017; %moment of inertia of internal mechanism
    I3 = .0006; %moment of inertia of the pendulum
    
    R = .15; %Radius of the outer spherical shell
    l = .12; %Lenght of the pendulum
    g = 9.81; %Gravity acceleration

    zeta = 0.03; %viscous damping coeficient associated with the pendulum-spher bearing
    
    gamma = 5*pi/180;%Plane slope (Simulation does not work with gamma = 0. Make gamma very small instead
    
    %auxiliar variables
    c1 = m3*R*l;
    c2 = m3*g*l;
    Mt = m1+m2+m3;
    
    %Diagonal members of M(q) //Notice they are constants
    m11 = Mt*R^2 + I1;
    m22 = m3*l^2 + I2 + I3;
    
    %State = x = [theta theta_dot phi phi_dot]
    
    %Initial State
    x0 = [0; 0; 0; 0];%initial state
    
    %Desired Initial and Final States for theta
    theta0d = [0;0];
    thetadf = [pi*2; 0];
    phid = asin(Mt*R*sin(gamma)/(m3*l));%equilibrium point
    
    
    %SIMULATION PARAMETERS
    tf =15; %Simulation end point (in seconds)
    nSimulationPoints = 500;
    
    %ode45 settings (for better results)
    options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4,1e-4, 1e-4, 1e-4, 1e-4,1e-4,1e-4]);
    
    %Extended State (Also store torque, xdesired and sliding mode deviation)
    %xx = [x0;theta_desired;theta_dot_desired;torque;S;s1;s2]
    xx0 = [x0;0;0;0;0;0;0];
    
    %Simulate the system without feedforward
    [T,X] = ode45(@(t,xx)sysODE(t,xx),linspace(0,tf,nSimulationPoints),xx0, options);
 
    figure();
    plot(T,X(:,1:6));
    hold on;
    plot(T,linspace(phid,phid,nSimulationPoints));
    hold off;
    legend('theta', 'theta dot', 'phi', 'phi dot','theta desired', 'theta dot desired', 'phi desired');
    xlabel('Time [s]');
    ylabel('Angle[Rad] and Angular Velocity[Rad/s]');
    %plot torque
    torquePlot = figure();
    torquePlot.Name = 'Torque Magnitude';

    %Hacking the ode to compute the torque, S, s1 and s2 values from their integrals
    UU = X(1:end,7:10)- circshift(X(1:end,7:10),[1,0]);
    TT = T(1:end) - circshift(T(1:end),[1,0]);
    
    U = UU(:,1)./TT;%Torque
    S = UU(:,2)./TT;
    s1 = UU(:,3)./TT;
    s2 = UU(:,4)./TT;

    %Plot Torque
    plot(T(2:end),U(2:end));
    ylabel('Torque [N.m]');
    xlabel('Time [s]');
    
    %Plot Sliding Mode Deviations
    figure();
    plot(T(2:end),S(2:end));
    hold on;
    plot(T(2:end),s1(2:end));
    plot(T(2:end),s2(2:end));
    hold off;
    legend('S', 's1', 's2');
    xlabel('Time [s]');
    ylabel('Sliding Variables Values');
    
    %x = [theta theta_dot phi phi_dot]
    function [xx_dot] = sysODE(t,xx)
        
        x = xx(1:4);%Extract state
        
        %Change state names to make it more readable
        theta = x(1);
        theta_dot = x(2);
        phi = x(3);
        phi_dot = x(4);
        
        %Compute desired trajectory
        thetad = planTrajectory(theta0d,thetadf,tf,t,1);%return [thetad,theta_dotd, theta_ddotd]
        xd = [thetad,phid,0,0];
        
        %M(q) components
        m12 = c1*cos(phi-gamma);
        m21 = m12;
        
        %Term used in N(q,q_dot)
        damp = zeta*(theta_dot+phi_dot);
        %N(q,q_dot) components
        n1 = damp + Mt*R*g*sin(gamma) - c1*sin(phi-gamma)*phi_dot^2;
        n2 = damp + c2*sin(phi);
        
        %Denominator of f(x) and b(x)
        f_den = m12*m21-m11*m22;
        b_den = -f_den;
        
        %Compute f(x) and b(x)
        f = [(m22*n1-m12*n2), (m11*n2-m21*n1)]./f_den;
        b = [m22-m12, m11-m21]./b_den;
        
        %Add disturbance
        d1 = rand()*1;
        d2 = rand()*1;
        
        neta = 20;%neta needs to be greater than alpha*d1_max + beta*d2_max
        [u,S,s1,s2] = computeSlidingModeTorque(x,xd,f,b,neta);
        
        %compute x_dot (theta_ddot and phi_ddot)
        x_dot(1) = x(2);
        x_dot(2) = f(1)+b(1)*u+d1;
        x_dot(3) = x(4);%phi_dot
        x_dot(4) = f(2)+b(2)*u+d2;
        
        xx_dot = [x_dot,thetad(2),thetad(3),u,S,s1,s2];
        
        xx_dot = xx_dot';

    end
      
simtimescale = 2; % by what factor to scale down the simulation playback wrt ODE solution

animeFigure = figure();
animeFigure.Name = 'Simulation';
i = 1;
j = 1;
mov(1:round(length(T)/simtimescale)) = struct('cdata', [],'colormap', []);
while i <=length(T)    %Display every frame
    %showHoopOnOtherHoop(animeFigure,[0,0],R1,X(i,1), X(i,5),lpend1,[0,0],R2,X(i,2), X(i,6),lpend2); % showSate plots the frame
    showHoop(animeFigure,[0,0],R,X(i,1), X(i,3),l,gamma); % showSate plots the frame
    %pause(1/(tf*100))
    % capture the frame in the animation in video object
    mov(j)=getframe(animeFigure);
    j= j+1;
    i= i+simtimescale;
end

%movie2avi(mov, 'Run2.avi', 'compression', 'none','quality',5)

end