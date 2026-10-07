function [funmin,minval,noreps] = SteepestDescentAlgorithm(x0,epsilon,gama,ogfun,varsvect)

    %x0 and varsvect defined as column vectors

    gradf = gradient(ogfun,varsvect);
    
    f_handle = matlabFunction(ogfun,'Vars',varsvect);
    gradf_handle = matlabFunction(gradf,'Vars',varsvect);
    
    funmin = [];
    funmin(:,end+1) = x0;
    funmininparg = num2cell(decimalAccuracyPoints(funmin(:,end),6));

    minval = [];
    minval(:,end+1) = feval(f_handle,funmininparg{:});
    noreps = 0;

    while norm(feval(gradf_handle,funmininparg{:}))>=epsilon
        noreps = noreps+1; 
        funmin(:,end+1) = funmin(:,end)-gama*feval(gradf_handle,funmininparg{:});
        funmininparg = num2cell(funmin(:,end));
        minval(:,end+1) = feval(f_handle,funmininparg{:});
    end
end