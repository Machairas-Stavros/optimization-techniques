function [funmin,minval,noreps] = SteepestDescentAlgorithmProjection(x0,epsilon,gama,sk,ogfun,varsvect,restrictions_boundaries)

    %x0 and varsvect defined as column vectors

    gradf = gradient(ogfun,varsvect);
    
    f_handle = matlabFunction(ogfun,'Vars',varsvect);
    gradf_handle = matlabFunction(gradf,'Vars',varsvect);

    restrictions = [varsvect>=restrictions_boundaries(:,1),varsvect<=restrictions_boundaries(:,2)];
    restrictions_handle = matlabFunction(restrictions,'Vars',varsvect);
    
    funmin = [];
    funmin(:,end+1) = x0;
    funmininparg = num2cell(decimalAccuracyPoints(funmin(:,end),6));
    
    minval = [];
    minval(:,end+1) = feval(f_handle,funmininparg{:});
    noreps = 0;

    while norm(feval(gradf_handle,funmininparg{:}))>epsilon
        noreps = noreps+1;

        xk_bar = funmin(:,end)-sk*feval(gradf_handle,funmininparg{:});
        rhinparg = num2cell(decimalAccuracyPoints(xk_bar,6));
        
        if ~all(feval(restrictions_handle,rhinparg{:}),'all')
            xk_bar = OrthogonalProjection(xk_bar,restrictions_boundaries,varsvect);
        end
            
        funmin(:,end+1) = funmin(:,end)+gama*(xk_bar-funmin(:,end));
        funmininparg = num2cell(decimalAccuracyPoints(funmin(:,end),6));
        minval(:,end+1) = feval(f_handle,funmininparg{:});
        
        if noreps>=500
            break
        end
    end
end