function xdot = f_est1(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: f_est1.m
%--------------------------------------------------------------------------
% Description: Flow map for discretized hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00
    
global s dimpsi dimt u A B
% state
retout = @(x,y) x(y);
psi = x(1:dimpsi);
t = x(end-2);
k = x(end-1);
j = x(end);

% differential equations
xdot = [A*psi + B*retout(u(t),size(u(t),end)); 1; 0; 0];
end