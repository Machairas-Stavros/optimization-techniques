function gauss_res = calculateGaussianLinearCombo(u1_grid,u2_grid,population)

    % This function calculates the approximation surface produced from each
    % chromosome.
    %
    % u1_grid and u2_grid is the result between the meshgrid implementation
    % between u1 and u2 vectors and represent the total number of
    % combinations between u1 and u2 points. Its size is N1xN2 where u1 has
    % N1 points and u2 N2 points.
    %
    % population has a size of NGx5xNSP
    %
    % ng is the number of gaussians used for the linear combination.
    %
    % nsp is the number of samples that exist in the population.
    %
    % gauss_res is a 3D matrix of size N1xN2xNSP including the estimated
    % surface for each chromosome.
    
    params = {};

    for i=1:size(population,2)
        params{end+1} = population(:,i,:);
    end

    gauss_res = createGaussianLinearCombo(u1_grid,u2_grid,params{1},params{2},params{3},params{4},params{5});

end