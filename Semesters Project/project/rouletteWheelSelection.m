function strongParents = rouletteWheelSelection(pop_mse,population)

    % This function implements the selection of chromosomes within the
    % population in a roulette wheel manner.
    % 
    % population is a 3D matrix of size NGx5xNSP
    %
    % ng is the number of gaussians used for the linear combination.
    %
    % nsp is the number of samples that exist in the population.
    %
    % pop_mse is the mean squared error for each chromosome within the
    % population. It has a size of either 1xK or Kx1 (doesn't really
    % matter).
    %
    % strongParents is the renewd population after the the selection
    % process.

    inverr = 1./pop_mse;

    inverr(isnan(inverr)) = 0;

    app_prop = inverr/sum(inverr);

    strongParents = population(:,:,randsample(size(population,3),size(population,3),true,app_prop));
end



