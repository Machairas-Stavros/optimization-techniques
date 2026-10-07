function minimum_interval = BisectDiffAlgorithm(search_interval,lamda,function123Diff)
    
    n = ceil((log10(search_interval(2)-search_interval(1))-log10(lamda))/log10(2));
    minimum_interval = search_interval;

    for k=1:n

        xk = (search_interval(1)+search_interval(2))/2;
        fxkdiff = function123Diff(xk);

        if fxkdiff==0
            search_interval(1) = xk;
            search_interval(2) = xk;
            minimum_interval(end+1,:) = search_interval;
            break;
        elseif fxkdiff>0
            search_interval(2) = xk;
        else
            search_interval(1) = xk;
        end

        minimum_interval(end+1,:) = search_interval;
    end

end