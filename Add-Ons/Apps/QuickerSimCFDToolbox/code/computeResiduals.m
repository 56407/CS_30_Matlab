function [stop, convergence] = computeResiduals(NS,F,u,sizep,convergence,maxres)

% computeResiduals - Compute velocity residuals and check stop criteria for
%                    a flow problem.
%
% This QuickerSim CFD Toolbox function is used to take current solution
% vector, right hand side vector and linear system matrix and compute
% current momentum equation residuals and check if the required convergence
% criteria have already been met.
% Mostly used to check convergence of the nonlinear term when iterating the
% nonlinearity in the Navier-Stokes equations.
%
% [stop, convergence] = computeResiduals(NS,F,u,sizep,convergence,maxres)
%
% Input arguments:
% NS - linearized matrix of the Navier-Stokes problem in the current
%      iteration. Matrix size must be n-by-n, where n equals total number
%      of unknowns in the mesh (2 or 3 velocity components times the number
%      of velocity nodes plus number of pressure nodes, i.e.:
%      dim*nVnodes+nPnodes, where dim = 2 or 3 denotes spatial dimension.
% F  - n-by-1 right hand side vector.
% u  - current solution vector.
% sizep - size (1-by-2 vector) of the p array storing coordinates of the
%      mesh nodes.
% convergence - at the beginning of the iteration process and at the first
%      call to computeResiduals function this can be an empty array [] and
%      during iterations with each call to computeResiduals function it is
%      suplemented with a new row which consists in subsequent columns of
%      the iteration number, residual of the x-velocity equation, residual
%      of the y-velocity equation, residual of the z-velocity equation and 
%      residual of the divergence equation. In 2D case the convergence
%      array does not include the column with the z-velocity residual.
% maxres - a scalar value specifying maximal residual value which is
%      acceptable for the iteration process to be broken. In other words: 
%      If all residuals drop below maxres, the function will set
%      output parameter stop equal to one which may be used outside the
%      function to stop the iteration process.
% 
% Output arguments:
% stop - a parameter set to 0 unless all residuals drop below maxres. In
%      the latter case the stop parameter is set to 1.
% convergence - a number_of_iterations-by-4 (in 2D case) or a
%      number_of_iterations-by-5 (in 3D case) array which stores history of
%      convergence. Each row corresponds to one iteration and columns store
%      iteration number, x-velocity residual, y-velocity residual,
%      z-velocity residual (only for 3D) and divergence of velocity
%      residual.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also PLOTRESIDUALS, INITSOLUTION.

nnodes = sizep(2);
dim = sizep(1);

res = F-NS*u;
resu = norm(res(1:nnodes));
resv = norm(res(nnodes+1:2*nnodes));

if(size(convergence,1)>0)
    iter = convergence(end,1)+1;
else
    iter = 0;
end

stop = 0;

if(dim == 2)
    resdivv = norm(res(2*nnodes+1:end));
    convergence = [convergence; iter resu resv resdivv];
    if(resu<maxres && resv<maxres && resdivv<maxres)
        stop = 1;
    end
elseif(dim == 3)
    resw = norm(res(2*nnodes+1:3*nnodes));
    resdivv = norm(res(3*nnodes+1:end));
    convergence = [convergence; iter resu resv resw resdivv];
    if(resu<maxres && resv<maxres && resw<maxres && resdivv<maxres)
        stop = 1;
    end
end



end