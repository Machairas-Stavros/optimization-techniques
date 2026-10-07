function population_mse = evaluatePopulation(u1_grid,u2_grid,population)

    % This function return the mse across each crhomosome of the
    % population.
    %
    % u1_grid and u2_grid is the grid of values containing all the possible
    % points combinations between the vectors u1 and u2. It derives from
    % the meshgrid function.
    %
    % population is the whole set of chromosomes. It is a 3D matrix of size
    % NGx5xNSP.
    %
    % NG is the total number of gaussian functions used.
    %
    % NSP is the total number of samples (chromosomes) inside the
    % population.
    % 
    % population_mse is the mse across each chromosome of the population.

    eval_funct = evaluationFunct(u1_grid,u2_grid);

    estim_funct = calculateGaussianLinearCombo(u1_grid,u2_grid,population);

    population_mse = mse3D(eval_funct,estim_funct);

    clear eval_funct estim_funct

end