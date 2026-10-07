function bin_pop = crossoverNmutation(population,nb,cb)

    % This function implements the crossover and mutation process.
    % Crossover is implemented between LSB and cb bits and mutation is
    % implemented randomly at only one bit each time placed between LSB and
    % cb.
    % 
    % population is a 3D matrix of size NGx5xNSP.
    %
    % NSP is the total number of samples that a population has.
    %
    % NG is the total number of gaussian functions that are used.
    %
    % cb represents the crossover bit, i.e., the bit up to which the
    % crossover process is implemented. Basically bits are considered to be
    % placed in a line from MSB to LSB where MSB is the (nb-1)th bit and
    % LSB is the 0th bit. In matlab MSB is the 1st element of the matrix 
    % and LSB is the last. So cb takes an integer value in the interval
    % {1,...,nb} and it is reverted.
    %
    % bin_pop is the population occured after the crossover and mutation
    % process.

    bin_pop = reshape(string(dec2bin(population)),size(population));

    bin_pop = char(bin_pop);

    for i=1:size(population,1)
        for j=1:size(population,2)
            bin_pop(i,nb-cb+1:size(bin_pop,2),j,1:size(population,3)) = bin_pop(i,nb-cb+1:size(bin_pop,2),j,randperm(size(population,3)));
        end
    end
    
    rmn = randsample(cb,1);
    rmn = nb-rmn+1;
    
    szbp = [size(bin_pop,1),1,size(bin_pop,3),size(bin_pop,4)];
    bin_pop(:,rmn,:,:) = reshape(dec2bin(bitxor(bin2dec(reshape(string(bin_pop(:,rmn,:,:)),szbp)),1)),szbp);
    bin_pop = string(bin_pop);
    bin_pop = bin2dec(bin_pop);
end