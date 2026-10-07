function minimum_interval = BisectAlgorithm(search_interval,epsilon,lamda,function123)

    minimum_interval = search_interval;

    while decimalAccuracyPoints(abs(search_interval(2)-search_interval(1)),10)>=lamda ...
            && decimalAccuracyPoints(abs(search_interval(2)-search_interval(1)),10)>2*epsilon

        x1 = (search_interval(2)+search_interval(1))/2-epsilon;
        x2 = (search_interval(2)+search_interval(1))/2+epsilon;
        fx1 = function123(x1);
        fx2 = function123(x2);

        if fx1<fx2
            search_interval(2) = x2;
        else
            search_interval(1) = x1;
        end
        minimum_interval(end+1,:) = search_interval;
    end

end