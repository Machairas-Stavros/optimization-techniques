function rounded_number = decimalAccuracyPoints(number,accuracy_points)
    
    %This returns a number keeping only the first accuracy_points decimal 
    % points of it.
    
    multiplier = 10^accuracy_points;
    rounded_number = round(number*multiplier)/multiplier;
end