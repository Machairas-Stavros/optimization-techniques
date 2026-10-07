function error = mse3D(x,y)

    % This function calculates the mse between two surfaces.
    % 
    % x is a 2D matrix of size NxM while y is a 3D matrix of size NxMxK.
    % error is the mean squared error across the 3rd dimension and is a 
    % matrix of size either 1xK or Kx1 (doesn't really matter).
    %
    % error represents the mean squared error between the surfaces.

    error = squeeze(sum((x-y).^2,[1 2])/(size(x,1)*size(x,2)));

end