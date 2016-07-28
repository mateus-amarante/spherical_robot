function [Xdot] = ode2Hoops(t,X)
global ms1 ms2 mp1 mp2 Is1 Is2 R1 R2 lpend1 lpend2 G zeta1 zeta2
% initialize the system variables
theta1 = X(1);
theta2= X(2);
theta1_dot = X(3);
theta2_dot = X(4);
phi1 = X(5);
phi2 =  X(6);
phi1_dot = X(7);
phi2_dot = X(8);

% define reference state variables
tstart=0;
tf = 10;
B = [0*pi/180;0;0;-360*pi/180;0;0];%thetaInit, thetaDotInit,thetaDdotInit,thetaFin, thetaDotFin,thetaDdotFin,
traj = getTrajecteoryPoint(t,tstart,tf,B);

theta2d =  traj(1) ;
theta1d = 0*pi/180-theta2d*R2/R1; 


theta1_dotd = -traj(2)*R2/R1;
theta2_dotd = traj(2);
gamma = -theta2 -theta1/R2 * R1;
phi1d = 0;%asin(sin(gamma)*((ms1+mp1)*R1/(mp1*lpend1)));
phi2d = 0;
phi1_dotd = 0;
phi2_dotd = 0;

% Get the inertia metrix
M_num = [
    [                                                                     Is1 + (mp1*((2*R1^2*cos(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2^2 + (2*R1^2*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2^2))/2 + (ms1*((2*R1^2*cos(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2^2 + (2*R1^2*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2^2))/2, Is1 + (mp1*((2*R1*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2 + (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)))/R2))/2 + (ms1*((2*R1*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2 + (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)))/R2))/2, (mp1*((2*R1*lpend1*cos(theta2 + (R1*theta1)/R2)*cos(phi1)*(R1 + R2))/R2 - (2*R1*lpend1*sin(theta2 + (R1*theta1)/R2)*sin(phi1)*(R1 + R2))/R2))/2,                                                         0]
    [ Is1 + (mp1*((2*R1*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2 + (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)))/R2))/2 + (ms1*((2*R1*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2)/R2 + (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)))/R2))/2,                                                                               Is1 + Is2 + (mp1*(2*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2 + 2*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2))^2))/2 + (ms1*(2*sin(theta2 + (R1*theta1)/R2)^2*(R1 + R2)^2 + 2*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2))^2))/2 + R2^2*mp2 + R2^2*ms2,          (mp1*(2*lpend1*cos(phi1)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) - 2*lpend1*sin(theta2 + (R1*theta1)/R2)*sin(phi1)*(R1 + R2)))/2,                                   R2*lpend2*mp2*cos(phi2)]
    [                                                                                                                                                                                         (mp1*((2*R1*lpend1*cos(theta2 + (R1*theta1)/R2)*cos(phi1)*(R1 + R2))/R2 - (2*R1*lpend1*sin(theta2 + (R1*theta1)/R2)*sin(phi1)*(R1 + R2))/R2))/2,                                                                                                                                                                                                  (mp1*(2*lpend1*cos(phi1)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) - 2*lpend1*sin(theta2 + (R1*theta1)/R2)*sin(phi1)*(R1 + R2)))/2,                                                                                       (mp1*(2*lpend1^2*cos(phi1)^2 + 2*lpend1^2*sin(phi1)^2))/2,                                                         0]
    [                                                                                                                                                                                                                                                                                                                                       0,                                                                                                                                                                                                                                                                                                                 R2*lpend2*mp2*cos(phi2),                                                                                                                                               0, (mp2*(2*lpend2^2*cos(phi2)^2 + 2*lpend2^2*sin(phi2)^2))/2]
    ];
CmT_num = [
    
- zeta1*(phi1_dot + theta1_dot) - theta2_dot*((mp1*((2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 - (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 - (2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 + (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2))/2 + (ms1*((2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 - (2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2))/2) + (ms1*(2*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 - theta1_dot*((ms1*((2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 + (2*R1^2*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2^2 - (2*R1^2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2^2 - (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2))/2 + (mp1*((2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 - (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 - (2*R1^2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2^2 + (2*R1^2*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2^2))/2) + (mp1*(2*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + (mp1*phi1_dot*((2*R1*lpend1*phi1_dot*cos(theta2 + (R1*theta1)/R2)*sin(phi1)*(R1 + R2))/R2 + (2*R1*lpend1*phi1_dot*sin(theta2 + (R1*theta1)/R2)*cos(phi1)*(R1 + R2))/R2))/2 - (G*R1*mp1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2 - (G*R1*ms1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2
- zeta2*(phi2_dot + theta2_dot) + (mp1*(2*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + theta1_dot*((mp1*(2*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2))*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + (2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 - (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2))/2 + (ms1*(2*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2))*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - (2*R1*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2 + (2*R1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2))/R2))/2) + theta2_dot*((ms1*(2*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2))*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + 2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + (mp1*(2*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2))*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + 2*sin(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*cos(theta2 + (R1*theta1)/R2)*(R1 + R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2) + (ms1*(2*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) - 2*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + (mp1*phi1_dot*(2*lpend1*phi1_dot*sin(phi1)*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + 2*lpend1*phi1_dot*sin(theta2 + (R1*theta1)/R2)*cos(phi1)*(R1 + R2)))/2 - G*mp1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - G*ms1*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + R2*lpend2*mp2*phi2_dot^2*sin(phi2)
- (mp1*(2*lpend1*phi1_dot*sin(phi1)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + 2*lpend1*phi1_dot*cos(phi1)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 - zeta1*(phi1_dot + theta1_dot) + (mp1*theta1_dot*(2*lpend1*sin(phi1)*((R1^2*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + 2*lpend1*cos(phi1)*((R1^2*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2^2 + (R1*theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + (mp1*theta2_dot*(2*lpend1*sin(phi1)*(theta2_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + 2*lpend1*cos(phi1)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + (mp1*phi1_dot*(2*lpend1*cos(phi1)*(theta2_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2) - lpend1*phi1_dot*sin(phi1) + (R1*theta1_dot*sin(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2) + 2*lpend1*sin(phi1)*(theta2_dot*(R2 + cos(theta2 + (R1*theta1)/R2)*(R1 + R2)) + lpend1*phi1_dot*cos(phi1) + (R1*theta1_dot*cos(theta2 + (R1*theta1)/R2)*(R1 + R2))/R2)))/2 + G*lpend1*mp1*sin(phi1)
- zeta2*(phi2_dot + theta2_dot) + (mp2*(2*lpend2^2*phi2_dot^2*cos(phi2)*sin(phi2) - 2*lpend2*phi2_dot*sin(phi2)*(R2*theta2_dot + lpend2*phi2_dot*cos(phi2))))/2 + (mp2*phi2_dot*(2*lpend2*sin(phi2)*(R2*theta2_dot + lpend2*phi2_dot*cos(phi2)) - 2*lpend2^2*phi2_dot*cos(phi2)*sin(phi2)))/2 + G*lpend2*mp2*sin(phi2)
];

%Position and velocity gains of Shell
kptheta1 = 11.5;
kvtheta1 = 11.5;
kptheta2 = 11.5;
kvtheta2 = 11.5;

%Position and velocity gains of Pendulum
kpphi1 = 30;
kvphi1 = 11;
kpphi2 = 30;
kvphi2 = 11;

% Controller - Computed torque control error matrix
e_ddot_1 = -kptheta1*(theta1-theta1d)  + -kvtheta1*(theta1_dot-theta1_dotd);
e_ddot_2 = -kptheta2*(theta2-theta2d)  + -kvtheta2*(theta2_dot-theta2_dotd);
e_ddot_3 = -kpphi1*(phi1-phi1d) + -kvphi1*(phi1_dot-phi1_dotd);
e_ddot_4 = -kpphi2*(phi2-phi2d) + -kvphi2*(phi2_dot-phi2_dotd);

e_ddot = [e_ddot_1; e_ddot_2; e_ddot_3; e_ddot_4];

tauVect = M_num*e_ddot - CmT_num;

% tau =[tauVect(1)-tauVect(3);tauVect(2)-tauVect(4);
%     tauVect(1)-tauVect(3);tauVect(2)-tauVect(4)];
tau =[tauVect(1);tauVect(2);
    tauVect(1);tauVect(2)];
%tau = [0; 0; 0; 0];
%Q_ddot =[theta1_ddot;theta2_ddot;phi1_ddot;phi2_ddot] ;
Q_ddot = M_num\(CmT_num+tau);
Xdot = [X(3);X(4);Q_ddot(1);Q_ddot(2);X(7);X(8);Q_ddot(3);Q_ddot(4)];


    function [Xd] = getTrajecteoryPoint(t,tstart,tf,B)
        %tstart = 0;
        tend = tf*1;
        %keep computing trajectory values for some time, stop few seconds
        %before the simulation ends
      
            %quintic polynimial profile
            A = [1 tstart tstart^2 tstart^3 tstart^4 tstart^5
                0 1 2*tstart 3*tstart^2 4*tstart^3 5*tstart^4
                0 0 2 6*tstart 12*tstart^2 20*tstart^3
                1 tend tend^2 tend^3 tend^4 tend^5
                0 1 2*tend 3*tend^2 4*tend^3 5*tend^4
                0 0 2 6*tend 12*tend^2 20*tend^3];
            %B = [90*pi/180;0;0;180*pi/180;0;0];
            
            const = inv(A)*B;
            q = const(1) + const(2)*t + const(3)*t^2 + const(4)*t^3 +const(5)*t^4 + const(6)*t^5;
            q_dot = const(2) + 2*const(3)*t + 3*const(4)*t^2 + 4*const(5)*t^3 + 5*const(6)*t^4;
            q_ddot =  2*const(3) + 6*const(4)*t + 12*const(5)*t^2 + 20*const(6)*t^3;
            
            %Xd = [theta, theta_dot, theta_ddot, phi, phi_dot, phi_ddot]
            Xd = [q, q_dot, q_ddot, 0, 0, 0];
            
       
        
    end
% X = theta1,theta2,theta1_dot,theta2_dot,phi1,phi2,phi1_dot,phi2_dot
end