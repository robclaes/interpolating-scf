function [v,l,hist] = SCF_monomial(A, v0, options)
    % SCF_monomial   Method to find smallest eigenvalue and corresponding
    % eigenvector of NEPv A(v). Based on standard SCF, but a matrix
    % polynomial in monomial basis is build using a number of previous points.
    %
    % [v,l,hist] = SCF_monomial(A,v0 [,options])
    %
    % options(struct) are
    %   max_iter: maximum number of SCF iterations
    %   tol: tolerance, required maximum norm of residual
    %   p_hist: maximum number of previous point to use for interpolation
    %   iter_add: number of smallest eigenpairs added in each iteration
    %   shift: select eigenvalues closest to shift. Ignored if NaN.
    % returns
    %   v: the eigenvector of the last SCF iterations
    %      (with smallest residual, if iter_add>1)
    %   l: the corresponding eigenvalue
    %   hist: struct containing the residual history
    max_iter = 100;
    tol = 1e-14;
    p_hist = 2;
    iter_add = 1;
    shift = NaN;
    n = length(v0);
    AA = @(l,v) A(v) - l*eye(n);
    prev_shift = false;
    if nargin == 3
        if isfield(options, 'max_iter')
            max_iter = options.max_iter;
        end
        if isfield(options, 'p_hist')
            p_hist = options.p_hist;
        end
        if isfield(options, 'iter_add')
            if (options.iter_add<=p_hist && rem(p_hist,iter_add)==0)
                iter_add = options.iter_add;
            else
                error('Wrong argument value in options: iter_add');
            end
        end
        if isfield(options, 'shift')
            shift = options.shift;
        end
        if isfield(options, 'tol')
            tol = options.tol;
        end
        if isfield(options, 'prev_shift')
            if (options.prev_shift && ~isnan(shift))
                error('prev_shift cannot be true if a shift is given')
            else
                prev_shift = options.prev_shift;
            end
        end
    end
    
    k = p_hist/iter_add;
    vhist = v0*ones(1,iter_add);
    lhist = [];
    res = [];
    for i = 1:k
        v = vhist(:,end-iter_add+1);
        [vv,ll] = eigs(A(v),iter_add, 'smallestreal');
        vv  = vv ./ vecnorm(vv);
        vhist = [vhist vv];
        lhist = [lhist diag(ll)];
        res = [res get_residuals(A,vv,diag(ll))];
    end
    
    i = k+1;
    while (max(res(:,end))> tol && size(res,2)<max_iter)
        vv = vhist(:,end-p_hist+1:end);
        ll = lhist(end-p_hist+1:end);
        if min(abs(diff(ll)))<1e-14
            break
        end
        [C0,C1] = monomial_polynomial(AA,ll,vv);
        [vv,ll] = smallest_eig(C0,C1,iter_add,shift,prev_shift,ll);
        vv = find_eigenvectors(vv,n);
        vhist = [vhist vv];
        lhist = [lhist ll];
        res = [res get_residuals(A,vv,ll)];
        i = i+1;
    end
    
    hist = struct;
    hist.res = res;
    hist.l = lhist;
    hist.v = vhist;
    
    [~,ind] = min(lhist(:,end));
    v = vv(:,ind);
    l = ll(ind);
end

function [v,l] = smallest_eig(A,B,k,shift,prev_shift,ll)
    % Find the smallest k eigenvalues of (A,B) and the corresponding
    % eigenvectors
    As = sparse(A);
    Bs = sparse(B);
    nbeigs = min(size(A,1),k);
    [V,D] = eigs(As,Bs,nbeigs,ll(end));
    L = diag(D);
    L_filter = logical(abs(imag(L))<1e-8);
    V = V(:,L_filter);
    L = real(L(L_filter));
    v = V(:,1);
    l = real(L(1));
end

function res =  get_residuals(A,vv,ll) 
    % GET_RESIDUALS(A,vv,ll) returns the residual for each pair 
    % (vv(:,i),ll(i)) in a column vector RES.
    res = zeros(size(ll));
    vv = vv ./ vecnorm(vv);
    for i = 1:length(ll)
        res(i) = norm( A(vv(:,i))*vv(:,i) - ll(i)*vv(:,i));
    end
end

function V = find_eigenvectors(vv,n)
    V = zeros(n,size(vv,2));
    for i = 1:size(vv,2)
        vi = reshape(vv(:,i),n,[]);
        [vvv,~,~] = svd(vi);
        V(:,i) = vvv(:,1);
    end
end