function [z, info] = BulirschStoer(dynFun,t,z0,tol)
% [z, info] = BulirschStoer(dynFun,t,z0,tol)
%
% Solves an initial value problem using the Bulirsch-Stoer method. This
% method is ideal for high-accuracy solutions to smooth initial value
% problems. 
% 
% Computes z(t) such that dz/dt = dynFun(t,z), starting from the initial
% state z0. The solution at the grid-points will be accurate to within tol.
%
% If the provided grid is insufficient, this function will automatically
% introduce intermediate grid points to achieve the required accuracy.
%
% INPUTS:
%   dynFun = function handle for the system dynamics
%       dz = dynFun(t,z)
%           t = scalar time
%           z = [nz,1] = state as column vector
%           dz = [nz,1] = derivative of state as column vector
%   t = [1,nt] = time grid-point vector
%   z0 = [nz,1] = initial state vector
%   tol = [nz,1] = error tolerance along each dimension. If tol is a
%       scalar, then all dimensions will satisfy that error tolerance.
%
% OUTPUTS:     (nt = n+1)
%   z = [nz,nt] = solution to the initial value problem
%
% NOTES:
%   Implementation details:
%   http://web.mit.edu/ehliu/Public/Spring2006/18.304/implementation_bulirsch_stoer.pdf
%


nt = length(t);
nz = size(z0,1);
z = zeros(nz,nt);
z(:,1) = z0;

info.error = zeros(size(z));
info.nFunEval = zeros(1,nt);
for i=2:nt
    tSpan = [t(i-1), t(i)];
    [zF, stepInfo] = BulirschStoerStep(dynFun, tSpan, z(:,i-1), tol);
    if strcmp(stepInfo.exit,'converged')  %Successful step!
       z(:,i) = zF; 
       info.error(:,i) = stepInfo.error;
       info.nFunEval(i) = stepInfo.nFunEval;
    else  %Failed to converge -- try again on a better mesh
        nGridRefine = 4;  % 
        time = linspace(tSpan(1),tSpan(2),nGridRefine);
        [zTmp, infoTmp] = BulirschStoer(dynFun,time,z0,tol);
        z(:,i) = zTmp(:,end);
        info.error(:,i) = infoTmp.error(:,end);
        info.nFunEval(i) = sum(infoTmp.nFunEval);
    end
end








end