function [oneHotDir, nfree] = g(x)
dims = size(x,2);
% Check all directions
oneHotDirs = [eye(dims); -eye(dims)];
possibleDirs = zeros(2*dims, dims);
for i = 1:2*dims
    % If direction is free, add direction to possibleDirs
    if isempty(intersect(x, x(end,:) + oneHotDirs(i,:), 'rows'))
        possibleDirs(i,:) = oneHotDirs(i,:);
    end
end
possibleDirs = possibleDirs(~all(possibleDirs == 0,2), :);

% Count how many of them are free
nfree = size(possibleDirs,1);
% Choose a free direction at random

if nfree == 0
   oneHotDir = zeros(1,dims); 
elseif nfree == 1
   oneHotDir = possibleDirs;
else
   oneHotDir = possibleDirs(randi(nfree), :);
end
end

