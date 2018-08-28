function [phi, convergence] = initScalarSolution(p,scalarValue)

% initScalarSolution - Set initial guess to the scalar field.
%
% This QuickerSim CFD Toolbox function initializes scalar solution field
% with a constant value of an initial guess and initializes convergence
% array to be used in computeScalarResiduals function.
%
% [phi, convergence] = initScalarSolution(p,scalarValue)
%
% Input arguments:
% p - array storing coordinates of all nodes in the mesh (detailed
%     description in help for importMeshGmsh function).
% scalarValue - a scalar with initial guess to the scalar solution field.
%
% Output arguments:
% phi - initialized solution vector.
% convergence - empty array ready to be passed to the computeScalarResidual
%       function.
% 
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: COMPUTESCALARRESIDUALS, PLOTSCALARRESIDUALS.

phi = scalarValue*ones(size(p,2),1);
convergence = [];

end