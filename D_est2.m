function inside = D_est2(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: D_est2.m
%--------------------------------------------------------------------------
% Description: Jump set for DiscHyGradAlg
% Return 0 if outside of D, and 1 if inside D
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

global kj1
% state
t = x(end-2);

if (kj1(t,1) ~= kj1(t+1,1) || kj1(t,1) == kj1(t+1,1)) && t < size(kj1,1)-1
    inside = 1;
else
    inside = 0;
end
end