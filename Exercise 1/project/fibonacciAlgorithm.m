function minimum_interval = fibonacciAlgorithm(search_interval,epsilon,lamda,function123)

    fn_lowerlimit = abs(search_interval(2)-search_interval(1))/lamda;
    find_n = false;
    n = 0;

    while ~find_n

        n = n+1;
        
        if fibonacci(n)>fn_lowerlimit
            find_n = true;
        end
    end

    minimum_interval = search_interval;
    
    x1 = search_interval(1)+(fibonacci(n-2)/fibonacci(n))*(search_interval(2)-search_interval(1));
    x2 = search_interval(1)+(fibonacci(n-1)/fibonacci(n))*(search_interval(2)-search_interval(1));

    for k=1:n

        fx1 = function123(x1);
        fx2 = function123(x2);

        if fx1>fx2
            search_interval(1) = x1;
            x1 = x2;
            x2 = search_interval(1)+(fibonacci(n-k-1)/fibonacci(n-k))*(search_interval(2)-search_interval(1));
        else
            search_interval(2) = x2;
            x2 = x1;
            x1 = search_interval(1)+(fibonacci(n-k-2)/fibonacci(n-k))*(search_interval(2)-search_interval(1));
        end
        
        minimum_interval(end+1,:) = search_interval;

        if k==n-2
            break
        end
    end

    x1n = x1;
    x2n = x1+epsilon;
    if function123(x1n)>function123(x2n)
        search_interval(1) = x1n;
    else
        search_interval(2) = x2n;
    end

    minimum_interval(end+1,:) = search_interval;

end


