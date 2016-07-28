%% Trajectory planning using 3rd order polynomial functions.
function [qd] = planTrajectory(q0,qf,tf,t,nofigure)

if(t>tf*.7)
    qd = [qf(1),0,0];
    return
end

tf = tf*.7;

M= [1 0 0 0;
    0 1 0 0;
    1 tf tf^2 tf^3;
    0 1 2*tf 3*tf^2];
b=[q0;qf];
a=M\b;

q = a(1)+a(2)*t+a(3)*t^2+a(4)*t^3;
q_dot  = a(2)+2*a(3)*t+3*a(4)*t^2;
q_ddot = 2*a(3)+6*a(4)*t;

qd = [q,q_dot,q_ddot];
    
if nofigure==1
    return
else

t=0:0.01:tf;    
    
figure('Name','Position (degree)');
plot(t,a(1)+a(2)*t+ a(3)*t.^2+a(4)*t.^3,'LineWidth',3);
title('Position (degree)')
grid

figure('Name','Velocity (degree/s)');
plot(t,a(2)*t+ 2*a(3)*t +3*a(4)*t.^2,'LineWidth',3);
title('Velocity (degree/s)')
grid

figure('Name','Acceleration (degree/s^2)');
plot(t, 2*a(3) +6*a(4)*t,'LineWidth',3);
title('Acceleration (degree/s^2)')
grid


end