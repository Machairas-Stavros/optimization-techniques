function [stop_iter,stop_reas,chromosome,min_pop_mse] = endcheck(max_reps,rep,epsilon,pop_mse,population,min_pop_mse)

    % This functions checks whether it is time to end the whole iteration 
    % process or not.
    %
    % max_reps is a variable indicating the maximum number of iterations
    % that the algorithm is enabled to execute.
    %
    % rep is the number of current iteration.
    %
    % epsilon is the desired accuracy.
    %
    % pop_mse is a vector containing the mse of each chromosome of the
    % whole population.
    %
    % population is the whole population given in order to select the best
    % chromosome.
    %
    % min_pop_mse is the minimum mse achieved up to the time of execution
    % over the whole iteration process.
    %
    % stop_iter indicated wehter it is time to end the operation.
    %
    % stop_reas indicated the reason causing the end of the algorithm.
    % Either "E": accuracy achieved or "R": repatitions exceeded. It is
    % empty if the algorithm is not deemed to stop at this point.
    %
    % chromosome is the chromosome selected if it is time to end the
    % operation or minimum value of mse is achieved.
    %
    % min_pop_mse is the renewed min_pop_mse if better accuracy is
    % achieved.
    
    if rep>max_reps || min(pop_mse)<epsilon
        stop_iter = true;
    else
        stop_iter = false;
    end

    if stop_iter==true
        if min(pop_mse)<epsilon
            stop_reas = "E";
        else
            stop_reas = "R";
        end
    else
        stop_reas = [];
    end

    if stop_iter==true
        chromosome = population(:,:,find(pop_mse==min(pop_mse)));
        if min(pop_mse)<min_pop_mse
            min_pop_mse = min(pop_mse);
        end
    elseif stop_iter==false && min(pop_mse)<min_pop_mse
        min_pop_mse = min(pop_mse);
        chromosome = population(:,:,find(pop_mse==min_pop_mse));
    else
        chromosome = [];
    end
end