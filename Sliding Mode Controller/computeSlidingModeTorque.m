function [u,S,s1,s2] = computeSlidingModeTorque(x,xd,f,b,neta)
    lambda1 = 8;
    lambda2 = 25;    
    alpha = 12;
    beta = 2;
    k = .5;
    kk = 40;
    %neta = 15;

        theta = x(1); theta_dot = x(2);% theta_ddot = x(3);
        thetad = xd(1); theta_dotd = xd(2); theta_ddotd = xd(3);
        
        phi = x(3); phi_dot = x(4);% theta_ddot = x(3);
        phid = xd(4); phi_dotd = xd(5); phi_ddotd = xd(6);
        
        e1 = theta-thetad;
        ed1 = theta_dot-theta_dotd;
        e2 = phi-phid;
        ed2 = phi_dot-phi_dotd;
        
        s1 = lambda1*e1 + ed1;
        s2 = lambda2*e2 + ed2;
        
        S = alpha*s1 + beta*s2;
        
        ueq1 = -(f(1) + lambda1*ed1 - theta_ddotd)/b(1);
        ueq2 = -(f(2) + lambda2*ed2 - phi_ddotd)/b(2);

        usw = (-beta*b(2)*ueq1 -alpha*b(1)*ueq2-neta*tanh(kk*S)-k*S)/(alpha*b(1)+beta*b(2));
        
        u = ueq1+ueq2+usw;
end