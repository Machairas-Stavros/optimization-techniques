function pop_ind = values2indexes(population,m1,m2,s1,s2,c,dap)

    % This function replaces the real values of m1,m2,s1,s2,c with indexes
    % in order to help the tranfromation of decimal numbers to binary.
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
    % NGx5xNSP containing values for m1,m2,s1,s2,c.
    %
    % ng is the number of gaussians used.
    %
    % nsp is the number of samples inside the population.
    %
    % dap are the decimal accuracy points used for the mapping process
    %
    % pop_ind is the population given back with each real values replaced
    % with indexes.

    pop_ind = zeros(size(population));

    params = cat(2,m1,m2,s1,s2,c);

    for i=1:size(params,2)
        param_map = containers.Map(decimalAccuracyPoints(params(:,i),dap),0:numel(params(:,i))-1);
        pop_ind(:,i,:) = cellfun(@(x) param_map(x), num2cell(decimalAccuracyPoints(population(:,i,:),dap)));
    end

    clear params param_map
end

