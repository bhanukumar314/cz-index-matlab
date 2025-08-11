% xdot = crtbpstm_ham(t,xx,mu)
% 
% first order circular restricted three-body equations of motion for Hamiltonian formulation
%
% input:
%  t = time (not required for autonomous crtbp system)
%  xx = state vector
%  mu = mass ratio
%
% output:
%  xdot = integration vector (derivatives)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = crtbp_ham(t,xx,mu)

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

xdot = [xdot;
        ydot;
        zdot; 
        pxdot;
        pydot
        pzdot];

end