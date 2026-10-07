function population = createPopulation(m1,m2,s1,s2,c,ng,nsp)
    
    % This function initiallizes the population on a random manner taking 
    % into account all the possbile values for each parameter.
    % 
    % m1 is a set of all the possible mean values of u1 vector.
    %
    % m2 is a set of all the possible mean values of u2 vector.
    %
    % s1 is a set of all the possible sigma values of u1 vector.
    %
    % s2 is a set of all the possible sigma values of u2 vector.
    %
    % c is a set of all the possible constant values.
    %
    % m1,m2,s1,s2,c are row vectors
    %
    % ng is the number of gaussians used.
    %
    % nsp is the number of samples inside the population.
    %
    % population is a 3D matrix of size NGx5xNSP.

    all_combs = combvec(m1,m2,s1,s2,c);

    all_combs = all_combs(:,randperm(size(all_combs,2)));

    ns = ng*nsp;

    population = all_combs(:,randsample(size(all_combs,2),ns,true));

    clear all_combs

    population = population(:,randperm(size(population,2)));

    population = reshape(population,5,ng,nsp);

    population = permute(population,[2 1 3]);
end
