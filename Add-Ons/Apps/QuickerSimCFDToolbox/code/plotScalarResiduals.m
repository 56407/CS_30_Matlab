function [] = plotScalarResiduals(convergence, figureId)

% plotScalarResiduals - Display convergence plot and print residuals to the
%                       console.
%
% This QuickerSim CFD Toolbox function, if called in each iteration of the
% solver run, displays actual convergence plot in a chosen Matlab window
% and prints iteration number and current residuals to the console.
%
% plotScalarResiduals(convergence, figureId)
%
% Input arguments:
% convergence - array initialized with initScalarSolution function and
%               filled in with the computeScalarResiduals function.
% figureId    - id of the Matlab figure window in which the plot should be 
%               displayed.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: COMPUTESCALARRESIDUALS, INITSCALARSOLUTION.

figure(figureId)
semilogy(convergence(:,1),convergence(:,2))
title('Passive scalar residual convergence plot');
xlabel('Iteration')
ylabel('Residual value [-]');
grid on;
drawnow;

disp(convergence(end,:));


end