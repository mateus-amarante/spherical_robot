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
function [] = showHoop(figureObj,HoopOffset,radHoop,thetaHoop, thetaPend,lenPend,gamma)
% switch the context figure to the passed parameter object;
figure(figureObj);
clf;
hold on

% plot circle for the shell
center = [radHoop*thetaHoop*cos(gamma)+HoopOffset(1),radHoop+HoopOffset(2)+radHoop*thetaHoop*sin(gamma)];
theta = linspace(0,2*pi,50);
Xo = center(1);
Yo = center(2);
CircX=radHoop*cos(theta)+Xo;
CircY=radHoop*sin(theta)+Yo;
plot(CircX,CircY,'color','r'); 

%Plot four points around circle that rotate with it to show rotational
%motion:
theta = linspace(0,2*pi,5);
CircX=radHoop*cos(theta-thetaHoop)+Xo;
CircY=radHoop*sin(theta-thetaHoop)+Yo;
plot(CircX,CircY,'.','color','r'); 

% plot pendulum
pendP1 = center;
pendP2 = center + [lenPend*sin(thetaPend),-lenPend*cos(thetaPend)];
lineArr = [pendP1;pendP2];
plot(lineArr(:,1),lineArr(:,2),'-o');

%plot the ground line
plot([0;radHoop*2*thetaHoop*cos(gamma)],[0;radHoop*2*thetaHoop*sin(gamma)])
%plot the max pendulum position
%plot([center(1);radHoop],[center(2);radHoop*2*thetaHoop*sin(gamma)])
axis equal
%axis([(floor(center(1)/(5*radHoop))-1)*5*radHoop; (floor(center(1)/(5*radHoop))+2)*5*radHoop; 0; 2*radHoop]);
axis([(floor(center(1)/(2*pi*radHoop)))*2*pi*radHoop-radHoop;(floor(center(1)/(2*pi*radHoop))+1)*2*pi*radHoop+radHoop; (floor(center(2)/(2*pi*radHoop)))*2*pi*radHoop-radHoop;(floor(center(2)/(2*pi*radHoop))+1)*2*pi*radHoop+radHoop;]);
hold off

end