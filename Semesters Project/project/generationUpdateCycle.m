function [population,pop_mse,endfunc,chromosome,min_pop_mse,stpr] = generationUpdateCycle(population,m1,m2,s1,s2,c,u1_grid,u2_grid,max_reps,k,e,min_pop_mse,nb,cb,dap)

    % This function comprises the whole evaluation, selection and evolution
    % process of the algorithm.
    % 
    % population is the whole set of chromosomes. It is a 3D matrix of size
    % NGx5xNSP.
    %
    % NG is the total number of gaussian functions used.
    %
    % NSP is the total number of samples (chromosomes) inside the
    % population.
    %
    % m1 is an 1xK matrix including the mean values across u1.
    % m2 is an 1xK matrix including the mean values across u2.
    %
    % s1 is an 1xK matrix including the sigma values across u1.
    % s2 is an 1xK matrix including the sigma values across u2.
    %
    % c is an 1xK matrix including the constant values for each
    % gaussian.
    %
    % m1,m2,s1,s2,c are column vectors where K is the total number of 
    % possible values that each vector can take.
    %
    % u1_grid and u2_grid is the result between the meshgrid implementation
    % between u1 and u2 vectors and represents the total number of
    % combinations between u1 and u2 points. Its size is N1xN2 where u1 has
    % N1 points and u2 N2 points.
    %
    % max_reps is the maximum number of repetitions up to its the algortihm
    % is permitted to run.
    %
    % k is the number of the current iteration.
    %
    % e is an upper limit for the mse. if mse reaches values belowe the 
    % predefined e then tha algorithm is allowed to be terminated.
    %
    % min_pop_mse is the minimum mse achieved up to this point.
    %
    % nb is the total number of bits used for the interpretation of
    % m1,m2,s1,s2,c.
    %
    % cb is represents the crossover bit.
    %
    % dap is a parameter used in the values2indexes. It defines the decimal
    % accuracy points for the mapping process.
    %
    % population is the renwed population after the crossover and mutation
    % process.
    %
    % pop_mse is the mse of the renewd population.
    %
    % endfunc indicated whether it is time to terminate the iteration
    % process.
    %
    % chromosome is the one combination of parameters achieving the minimum
    % mse.
    %
    % min_pop_mse is the renewd minimum mse achieved.
    %
    % stpr is the terminating reason of the alogrithm.

    if nargin==14
        dap = 5;
    end

    pop_mse = evaluatePopulation(u1_grid,u2_grid,population);

    [endfunc,stpr,chromosome,min_pop_mse] = endcheck(max_reps,k,e,pop_mse,population,min_pop_mse);
    
    if endfunc
        return
    end
    
    population = rouletteWheelSelection(pop_mse,population);
    
    population = values2indexes(population,m1,m2,s1,s2,c,dap);
    
    population = crossoverNmutation(population,nb,cb);
    
    population = indexes2values(population,m1,m2,s1,s2,c);
    
end