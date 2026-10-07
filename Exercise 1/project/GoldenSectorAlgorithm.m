function minimum_interval = GoldenSectorAlgorithm(search_interval,lamda,function123)

    gama = (sqrt(5)-1)/2;
    x1 = search_interval(1)+(1-gama)*(search_interval(2)-search_interval(1));
    x2 = search_interval(1)+gama*(search_interval(2)-search_interval(1));
    minimum_interval = search_interval;

    while decimalAccuracyPoints(abs(search_interval(2)-search_interval(1)),10)>=lamda

        fx1 = function123(x1);
        fx2 = function123(x2);

        if fx1>fx2
            search_interval(1) = x1;
            x1 = x2;
            x2 = search_interval(1)+gama*(search_interval(2)-search_interval(1));
        else
            search_interval(2) = x2;
            x2 = x1;
            x1 = search_interval(1)+(1-gama)*(search_interval(2)-search_interval(1));
        end
        minimum_interval(end+1,:) = search_interval;
    end

end