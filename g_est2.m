function xplus = g_est2(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: g_est2.m
%--------------------------------------------------------------------------
% Description: Jump map for DiscHyGradAlg
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

global s gamma_c1 gamma_d1 dimpsi dimth u H J D error theta psi kj1
% state
retout = @(x,y) x(y);
thetarows = x(1:end-3);
t = x(end-2); k = x(end-1); j = x(end);

for i = 1:dimth(1)
    theta2(i,:) = thetarows(1+(i-1)*dimth(2):dimth(2)+(i-1)*dimth(2));
end
psitrue1 = psi(t,1:dimpsi)'+psi(t,2+dimpsi:2+dimpsi*2-1)'; %2x1
psitrue = psi(t,:)';
error = (theta2-theta)*psitrue;

if kj1(t,1) ~= kj1(t+1,1)
    theta2 = theta2-s*gamma_c1*(psitrue1*psitrue1')*theta2;
end
if kj1(t,1) == kj1(t+1,1)
    psitrue2 = H*psitrue1+J*retout(u(t),size(u(t),end));
    theta2 = theta2-gamma_d1*psitrue2*psitrue2'*theta2/(1+gamma_d1*psitrue2'*psitrue2);
end

for i = 1:dimth(1)
    theta2rows(1+(i-1)*dimth(2):dimth(2)+(i-1)*dimth(2)) = theta2(i,:);
end

if kj1(t,1) ~= kj1(t+1,1)
xplus = [theta2rows'; t+1; k+s; j];
end
if kj1(t,1) == kj1(t+1,1)
xplus = [theta2rows'; t+1; k; j+1];
end
end