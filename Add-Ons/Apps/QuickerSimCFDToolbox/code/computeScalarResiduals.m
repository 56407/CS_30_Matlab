function  [stop, convergence] = computeScalarResiduals(D,F,T,convergence,maxres)

% computeScalraResiduals - Compute residuals of a scalar transport problem.
%
% This QuickerSim CFD Toolbox function computes current residual of a
% scalar transport problem.
%
% [stop, convergence] = computeScalarResiduals(D,F,T,convergence,maxres)
%
% Input arguments:
% D - global matrix of the scalar transport problem.
% F - right hand side vector.
% T - current solution vector.
% convergence - array with iteration number (first column) at residual
%     value (second column) of the iteration process. At the beginning
%     empty array initialized by the initScalarSolution function.
% maxres - a scalar value specifying maximal residual value to be achieved
%     for a converged solution.
% 
% Output arguments:
% stop - a parameter set to 0 unless the residual drops below maxres. In
%      the latter case the stop parameter is set to 1.
% convergence - input array supplemented with values from the new
%      iteration.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also INITSCALARSOLUTION, PLOTSCALARRESIDUALS.

if(size(convergence,1)>0)
    iter = convergence(end,1)+1;
else
    iter = 0;
end

stop = 0;
res = norm(D*T-F);
convergence = [convergence; iter res];

if(res<maxres)
    stop = 1;
end

end