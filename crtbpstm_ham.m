% xdot = crtbpstm_ham(t,xx,mu)
% 
% first order circular restricted three-body equations of motion + state
% transition matrix for Hamiltonian formulation
%
% input:
%  t = time (not required for autonomous pcrtbp system)
%  xx = state vector
%  mu = mass ratio
%
% output:
%  xdot = integration vector (derivatives)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = crtbpstm_ham(t,xx,mu)

x = xx(1);
y = xx(2);
z = xx(3);
px = xx(4);
py = xx(5);
pz = xx(6);

r1 = ((x+mu)^2 + y^2 + z^2)^0.5;
r2 = ((x-1+mu)^2 + y^2 + z^2)^0.5;

xdot = px + y;
ydot = py - x;
zdot = pz;
pxdot = py - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
pydot = -px + (- (1-mu)/r1^3 - mu/r2^3)*y;
pzdot = (- (1-mu)/r1^3 - mu/r2^3)*z;

phi=reshape(xx(7:6*6+6),6,6);

A = [                                                                                                                                                                                                            0, 1,  0, 1, 0, 0;
                                                                                                                                                                                                                -1, 0,  0, 0, 1, 0;
                                                                                                                                                                                                                 0, 0,  0, 0, 0, 1;
 (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(mu + x - 1)*(mu + x - 1))/(((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(mu + x)*(mu + x)*(mu - 1))/(((mu + x)^2 + y^2 + z^2)^(5/2)),      (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),       (3*z*mu*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),  0, 1, 0;
  y*((3*mu*(mu + x - 1))/(((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(mu + x)*(mu - 1))/(((mu + x)^2 + y^2 + z^2)^(5/2))),       (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + y*((3*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2)),       y*(3*mu*z/(((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(mu - 1))/(((mu + x)^2 + y^2 + z^2)^(5/2))), -1, 0, 0;
  z*((3*mu*(mu + x - 1))/(((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(mu + x)*(mu - 1))/(((mu + x)^2 + y^2 + z^2)^(5/2))),       z*(3*mu*y/(((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(mu - 1))/(((mu + x)^2 + y^2 + z^2)^(5/2))),       (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + z*((3*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2)), 0, 0, 0];

phidot = A*phi;

xdot = [xdot;
        ydot;
        zdot; 
        pxdot;
        pydot;
        pzdot];

xdot(7:6*6+6)=reshape(phidot,6*6,1);
