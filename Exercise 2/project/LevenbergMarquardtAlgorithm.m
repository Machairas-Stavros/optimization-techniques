function [funmin,noreps,glist,armparams] = LevenbergMarquardtAlgorithm(x0,epsilon,ogfun,varsvect,gama,armijoparams)
    
    %varsvect is a column vector with the system variables and x0 a column
    %vector with the initial values of the varsvect(1:end-1) variables

    if nargin==5 && isstring(gama) && isequal(gama,"Armijo")
        armijoparams = "Default";
    elseif nargin==5
        armijoparams = [];
    end

    if isstring(armijoparams) && armijoparams=="Default"
        alpha = 0.1;
        beta = 0.5;
        sigma = 1;
        armparams = [alpha,beta,sigma];
    elseif isstring(armijoparams) && armijoparams=="Random"
        alpha = 10^-5+rand*(0.1-10^-5);
        beta = 0.1+rand*(0.5-0.1);
        sigma = rand;
        armparams = [alpha,beta,sigma];
    elseif isnumeric(armijoparams) && ~isempty(armijoparams)
        alpha = armijoparams(1);
        beta = armijoparams(2);
        sigma = armijoparams(3);
        armparams = [alpha,beta,sigma];
    end

    glist = [];

    gradf = gradient(ogfun,varsvect(1:end-1));
    hessianf = hessian(ogfun,varsvect(1:end-1));

    f_handle = matlabFunction(ogfun,'Vars',varsvect(1:end-1));
    gradf_handle = matlabFunction(gradf,'Vars',varsvect(1:end-1));
    hessianf_handle = matlabFunction(hessianf,'Vars',varsvect(1:end-1));
    
    funmin = x0;
    noreps = 0;
    
    funmininparg = num2cell(decimalAccuracyPoints(funmin,6));

    while norm(feval(gradf_handle,funmininparg{:}))>=epsilon

        mk_pd = ceil(max(abs(eig(feval(hessianf_handle,funmininparg{:})))));
        deltakf = hessianf+mk_pd*eye(length(varsvect)-1);
        deltakf_handle = matlabFunction(deltakf,'Vars',varsvect(1:end-1));

        if all(eig(feval(deltakf_handle,funmininparg{:}))>0) && det(feval(deltakf_handle,funmininparg{:}))~=0
            noreps = noreps+1;
            
            dkf = -inv(deltakf)*gradf;
            dkf_handle = matlabFunction(dkf,'Vars',varsvect(1:end-1));

            dk = feval(dkf_handle,funmininparg{:});

            if isnumeric(gama) && isequal(size(gama),[1,1]) && gama>0 && gama<=1
                gamak = gama;
            elseif isstring(gama) && isequal(gama,"Minimize")
                modifiedfinparg = num2cell(varsvect(1:end-1)+varsvect(end)*dkf);
                modifiedf = feval(f_handle,modifiedfinparg{:});

                modifiedf_gkdiff = diff(modifiedf,varsvect(end));

                modifiedf_gkdiff_subs = subs(modifiedf_gkdiff,varsvect(1:end-1),decimalAccuracyPoints(funmin,4));
                gamak = double(solve(modifiedf_gkdiff_subs==0,varsvect(end),'Real',true));
                gamak = gamak(gamak>=0);
    
                if length(gamak)>1
                    if sum(gamak<=1)>0
                        gamak = gamak(gamak<=1);
                        selind = randsample(length(gamak),1);
                        gamak = gamak(selind);
                    else
                        gamak = min(gamak);
                    end
                elseif isempty(gamak)
                    return
                end

                glist(end+1,:) = gamak;
            elseif isstring(gama) && isequal(gama,"Armijo")
                mk_arm = 0;
                while true
                    gamak = sigma*beta^mk_arm;
                    tempfunmininparg = num2cell(funmin+gamak*dk);
                    if feval(f_handle,funmininparg{:})-feval(f_handle,tempfunmininparg{:})>= -alpha*gamak* ...
                            transpose(dk)*feval(gradf_handle,funmininparg{:})
                        break
                    end
                    mk_arm = mk_arm+1;
                end
            end

            funmin = funmin+gamak*dk;
            funmininparg = num2cell(funmin);
        else
            return
        end

    end