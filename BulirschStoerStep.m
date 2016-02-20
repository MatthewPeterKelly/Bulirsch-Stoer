function [zF, info] = BulirschStoerStep(dynFun, tSpan, z0, tol)

nRefineMax = 8;   % How many attempts at refinement before giving up

%%%%%%%%%%%%%%%%%%%%%%%%%

n = 2*(1:nRefineMax);

nz = size(z0,1);

if length(tol)==1
    tol = tol*ones(size(z0));
end

T = zeros(nz,nRefineMax,nRefineMax);   %Extrapolation table
E = zeros(nz,nRefineMax);   %Error estimate table

info.exit = 'maxRefine';  %Assume that we fail to meet tolerance
for j=1:nRefineMax
    
    % Compute the current estimate of the solution:
    [~,z] = modifiedMidpointRule(dynFun, tSpan, z0, n(j));
    T(:,j,1) = z(:,end);
    
    if j>1
        
        % Compute the extrapolation table entries:
        for k=2:j
            num = T(:,j,k-1) - T(:,j-1,k-1);
            den = (n(j)/(n(j-k+1)))^2 - 1;
            T(:,j,k) = T(:,j,k-1) + num/den;
        end
        
        % Compute the error estimates:
        E(:,j) = abs(T(:,j,j-1) - T(:,j,j));
        
        % Check convergence:
        if all(E(:,j)<tol)
            info.exit = 'converged';
            break;
        end
    end
    
end

% Other useful things:
info.error = E(:,j);     %Error estimate
info.nFunEval = sum(n(1:j));    %number of function evaluations
info.nRefine = j;

% Return the estimate of the solution:
zF = T(:,j,j);

end