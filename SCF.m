function [v, l, hist] = SCF(A, v0, options)
    % SCF   Method to find smallest eigenvalue and corresponding
    % eigenvector of NEPv A(v).
    %
    % [v,l,hist] = SCF(A,v0 [,options])
    %
    % options(struct) are
    %   max_iter: maximum number of SCF iterations
    %   tol: tolerance, required maximum norm of residual
    % returns
    %   v: the eigenvector of the last SCF iterations
    %   l: the corresponding eigenvalue
    %   hist: struct containing the residual history
    
    max_iter = 100;
    tol = 1e-14;
    target = 'smallestreal';
    if nargin == 3
        if isfield(options, 'max_iter')
            max_iter = options.max_iter;
        end
        if isfield(options, 'tol')
            tol = options.tol;
        end
        if isfield(options, 'target')
            target = options.target;
        end
    end
    vv = [];
    [v,l] = eigs(A(v0),1, target);
    v = v / norm(v);
    res = norm(A(v)*v-l*v);
    ll = l;
    vv = v;
    while( norm(res(end))>tol && length(res)<max_iter )
        if length(ll)<2
            [v,l] = eigs(A(v),1, target);
        else
            [v,l] = eigs(A(v),1,ll(end));
        end
        v = v / norm(v);
        res = [res norm(A(v)*v-l*v)];
        ll = [ll l];
        vv = [vv v];
    end
    hist = struct('res', res, 'l', ll,'v',vv);
end