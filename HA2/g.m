function [oneHotDir, nbFree] = g(x)
dims = length(x(0,:));
% Check all directions
oneHotDirs = [eye(dims); -eye(dims)];
possibleDirs = zeros(2*dims, dims);
for i = 1:2*dims
    % If direction is free, add direction to possibleDirs
    if ~isempty(intersect(x, x(end,:) + oneHotDirs(i,:), 'rows'))
        possibleDirs(i,:) = oneHotDirs(i,:);
    end
end
possibleDirs = possibleDirs(~all(possibleDirs == 0,2), :);
% Count how many of them are free
nbFree = length(possibleDirs);
% Choose a free direction at random
oneHotDir = possibleDirs(randi(length(possibleDirs)), :);
end

