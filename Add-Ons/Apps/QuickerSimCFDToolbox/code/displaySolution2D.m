function [] = displaySolution2D(p, t, field, label, varargin)

% displaySolution2D - Display scalar field in 2-D.
%
% This QuickerSim CFD Toolbox function displays a 2-D scalar field.
%
% displaySolution2D(p, t, field, label);
% displaySolution2D(p, t, field, label, colorBarLimits);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% field - column vector with scalar field to be displayed.
% label - a string which is displayed as the title over the color map.
% colorBarLimits - a two element vector defining minimum and maximum value
%     represented by deep blue and hot red colors in the colormap. All
%     field values below and above these specified values will be colored
%     with the color corresponding to minimum and maximum values.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: DISPLAYMESH2D, IMPORTMESHGMSH.

nPnodes = max(max(t(1:3,:)));

trisurf(t(1:3,:)',p(1,1:nPnodes)',p(2,1:nPnodes)',zeros(nPnodes,1),field(1:nPnodes))
shading('interp')
colormap('jet')
colorbar
axis equal
view([0 90])
title(label)
xlabel('x')
ylabel('y')

if(~isempty(varargin))
    caxis(varargin{1})
end

end