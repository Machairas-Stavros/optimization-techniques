function pop_val = indexes2values(population,m1,m2,s1,s2,c)

    % This function replaces the indexes of each parameter with their real
    % values respsect to the m1,m2,s1,s2,c field of values.
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
    % m1,m2,s1,s2,c are column vectors.
    %
    % population is a set of chromosomes. It is a 3D matrix of size 
    % NGx5xNSP containing indexes for m1,m2,s1,s2,c.
    %
    % ng is the number of gaussians used.
    %
    % nsp is the number of samples inside the population.
    %
    % pop_val is the population filled with values.

    pop_val = zeros(size(population));

    params = cat(2,m1,m2,s1,s2,c);

    for i=1:size(params,2)
        pop_val(:,i,:) = reshape(params(population(:,i,:)+1,i),size(population(:,i,:)));
    end

    clear params param_map
end
