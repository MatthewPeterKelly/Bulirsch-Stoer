function [t,z] = modifiedMidpointRule(dynFun, tSpan, z0, n)
% [t,z] = modifiedMidpointRule(dynFun, tSpan, z0, n)
%
% Approximates the solution to the initial value problem by numerical
% integration with the modified mid-point rule. 
%
% INPUTS:
%   dynFun = function handle for the system dynamics
%       dz = dynFun(t,z)
%           t = scalar time
%           z = [nz,1] = state as column vector
%           dz = [nz,1] = derivative of state as column vector
%   tSpan = [1,2] = [t0,tF] = time span for simulation
%   z0 = [nz,1] = initial state vector
%   n = scalar integer number of steps.   (require:  n > 2)
%
% OUTPUTS:     (nt = n+1)
%   t = [1,nt] = time stamps for intermediate points
%   z = [nz,nt] = state estimate at final time
%
% NOTES:
%   Implementation details:
%   http://web.mit.edu/ehliu/Public/Spring2006/18.304/implementation_bulirsch_stoer.pdf
%

nt = n+1;
nz = size(z0,1);

t0 = tSpan(1);
tF = tSpan(2);
h = (tF-t0)/n;
t = linspace(tSpan(1), tSpan(2), nt);
z = zeros(nz,n);

% Initialize and then run the modified mid-point method.
z(:,1) = z0;
z(:,2) = z0 + h*dynFun(t0,z0);
for i=3:nt
   z(:,i) = z(:,i-2) + 2*h*dynFun(t(i-1),z(:,i-1)); 
end

% Refine the final point using the dynamics function at final point.
z(:,nt) = 0.5*(z(:,nt) + z(:,nt-1) + h*dynFun(t(nt),z(:,nt)));

end