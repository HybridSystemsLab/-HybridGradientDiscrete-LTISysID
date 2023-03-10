function xplus = g_est1(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: g_est1.m
%--------------------------------------------------------------------------
% Description: Jump map for discretized hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

global s dimpsi u H J D
% state
retout = @(x,y) x(y);
psi = x(1:dimpsi);
t = x(end-2);
k = x(end-1);
j = x(end);

if t >= k+s
xplus = [psi; t; k+s; j];
end
if D(psi,t)
xplus = [H*psi+J*retout(u(t),size(u(t),end)); t; k; j+1];
end
end