function [] = plotResiduals(convergence, figureId)

% plotResiduals - Display convergence plot and print residuals to console.
%
% This QuickerSim CFD Toolbox function, if called in each iteration of the
% solver run, displays actual convergence of the velocity residuals in a
% chosen Matlab window and prints iteration number and current residuals to
% console window.
%
% plotResiduals(convergence, figureId)
%
% Input arguments:
% convergence - array initialized with initSolution function and filled in
%               at least once with the computeResiduals function.
% figureId    - id of the Matlab figure window in which the plot should be 
%               displayed.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: COMPUTERESIDUALS, INITSOLUTION.

figure(figureId)

if(size(convergence,2)==4) % dim = 2
    semilogy(convergence(:,1),convergence(:,2:4))
    legend('x-velocity residual','y-velocity residual','Continuity residual');
    tmpconv = convergence;
    tmpconv(tmpconv<1e-12) = mean(mean(tmpconv));
    ylim([min(min(tmpconv(:,2:3))), max(max(tmpconv(:,2:3)))])
else % dim = 3
    semilogy(convergence(:,1),convergence(:,2:5))
    tmpconv = convergence;
    tmpconv(tmpconv<1e-12) = mean(mean(tmpconv));
    ylim([min(min(tmpconv(:,2:4))), max(max(tmpconv(:,2:4)))])
    legend('x-velocity residual','y-velocity residual','z-velocity residual','Continuity residual');
end

%semilogy(convergence(:,1),convergence(:,2:3))
title('Residuals convergence plot');
xlabel('Iteration')
ylabel('Residual value [-]');
%legend('x-velocity residual','y-velocity residual');

grid on;
drawnow;

disp(convergence(end,:));

end