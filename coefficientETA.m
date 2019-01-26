function ret = coefficientETA(x)
    load('derenzo.mat', 'centroids');
    load('derenzo.mat', 'radius');
    
    xmax = max(centroids(:, 1));
    xmin = min(centroids(:, 1));
    ymax = max(centroids(:, 2));
    ymin = min(centroids(:, 2));

    r = max(xmax - xmin, ymax - ymin) * 1.15;

    centroids(:, 1) = (centroids(:, 1) - xmin*1.2) / r;
    centroids(:, 2) = (centroids(:, 2) - ymin*1.2) / r;


    radius = radius / r;

    ret = 0.1*ones(1, size(x, 2));

    for i = 1:size(centroids, 1)
        ret = ret + 16*radius(i)* ((x(1,:) - centroids(i, 1)).^2 + (x(2,:) - centroids(i, 2)).^2 < 2*radius(i)^2);
    end
end