%{
The showHoop function generates a plot of the current state of the single
hoop on flat surface system.

Parameters:
HoopOffset: a cartesian offset to be applied to the center of the hoop
(used if sliding motion is required)
radHoop: the hoop radius
thetaHoop: the angular displacement of the hoop in radians
thetaPend: the angular displacement of the pendulum in radians
lenPend: the length of the pendulum
%}
function [] = showHoopOnOtherHoop(figureObj,t_HoopOffset,t_radHoop,t_thetaHoop, t_thetaPend,t_lenPend,b_HoopOffset,b_Radius,b_thetaHoop, b_thetaPend,b_lenPend)
% switch the context figure to the passed parameter object;
figure(figureObj);
clf;
hold on
gamma = -t_radHoop/b_Radius * t_thetaHoop;

%screenCenter = []
% Plot the base:
centerL = [b_thetaHoop*b_Radius,b_Radius];
theta = linspace(0,2*pi,100);
Xo = centerL(1);
Yo = centerL(2);
CircX=b_Radius*cos(theta)+Xo;
CircY=b_Radius*sin(theta)+Yo;
plot(CircX,CircY,'color','g');
%Plot four points around circle that rotate with it to show rotational
%motion:
theta = linspace(0,2*pi,5);
CircX=b_Radius*cos(theta-b_thetaHoop)+Xo;
CircY=b_Radius*sin(theta-b_thetaHoop)+Yo;
plot(CircX,CircY,'.','color','g'); 
% plot pendulum
pendP1 = centerL;
pendP2 = centerL + [b_lenPend*sin(b_thetaPend),-b_lenPend*cos(b_thetaPend)];
lineArr = [pendP1;pendP2];
plot(lineArr(:,1),lineArr(:,2),'-o');


% plot circle for the shell
%center = [radHoop*thetaHoop*cos(gamma)+HoopOffset(1),radHoop+HoopOffset(2)+radHoop*thetaHoop*sin(gamma)];
center = [centerL(1)+ -(t_radHoop +b_Radius)*sin(gamma),centerL(2)+ (t_radHoop +b_Radius)*cos(gamma)];
theta = linspace(0,2*pi,50);
Xo = center(1);
Yo = center(2);
CircX=t_radHoop*cos(theta)+Xo;
CircY=t_radHoop*sin(theta)+Yo;
plot(CircX,CircY,'color','r'); 

%Plot four points around circle that rotate with it to show rotational
%motion:
theta = linspace(0,2*pi,5);
CircX=t_radHoop*cos(theta-t_thetaHoop)+Xo;
CircY=t_radHoop*sin(theta-t_thetaHoop)+Yo;
plot(CircX,CircY,'.','color','r'); 

% plot pendulum
pendP1 = center;
pendP2 = center + [t_lenPend*sin(t_thetaPend),-t_lenPend*cos(t_thetaPend)];
lineArr = [pendP1;pendP2];
plot(lineArr(:,1),lineArr(:,2),'-o');

%plot the base circle
%plot([0;radHoop*2*thetaHoop*cos(gamma)],[0;radHoop*2*thetaHoop*sin(gamma)])
%plot the max pendulum position
%plot([center(1);radHoop],[center(2);radHoop*2*thetaHoop*sin(gamma)])
axis equal
%axis([(floor(center(1)/(5*radHoop))-1)*5*radHoop; (floor(center(1)/(5*radHoop))+2)*5*radHoop; 0; 2*radHoop]);
axis([-1+(floor(center(1)/(2*pi*t_radHoop)))*2*pi*t_radHoop-t_radHoop;(floor(center(1)/(2*pi*t_radHoop))+1)*2*pi*t_radHoop+t_radHoop; (floor(center(2)/(2*pi*t_radHoop)))*2*pi*t_radHoop-t_radHoop;(floor(center(2)/(2*pi*t_radHoop))+1)*2*pi*t_radHoop+t_radHoop;]);
hold off

end