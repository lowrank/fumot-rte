function [ ret ] = coefficientXF( x )
%ABSORPTIONFF Absorption coef for fluorescence material.
%   x coordinate, vectorized.
%   ret ~0.2
N = 500;
[X,Y] =meshgrid(linspace(-0.5,0.5, N));
X = X(:);
Y = Y(:);
ph = phantom('Modified Shepp-Logan', N);
F = scatteredInterpolant(X,Y,ph(:));
ret = 0.1 * F(x(1,:)-0.5, x(2,:)-0.50) + 0.05;
% ret = 0.2 + 0.05 * x(1,:);

% load derenzo.mat
% 
% xmax = max(centroids(:, 1));
% xmin = min(centroids(:, 1));
% ymax = max(centroids(:, 2));
% ymin = min(centroids(:, 2));
% 
% r = max(xmax - xmin, ymax - ymin) * 1.15;
% 
% centroids(:, 1) = (centroids(:, 1) - xmin*1.2) / r;
% centroids(:, 2) = (centroids(:, 2) - ymin*1.2) / r;
% 
% 
% radius = radius / r;
% 
% ret = zeros(1, size(x, 2));
% 
% for i = 1:size(centroids, 1)
%     ret = ret + 80*radius(i)* ((x(1,:) - centroids(i, 1)).^2 + (x(2,:) - centroids(i, 2)).^2 < 2*radius(i)^2);
% end


end
