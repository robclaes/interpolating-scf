function [v, l, hist] = SCF_DIIS(A, v0, options)
    % SCF_DIIS   Method to find smallest eigenvalue and corresponding
    % eigenvector of NEPv A(v). SCF iterations with application of DIIS
    %
    % [v,l,hist] = SCF(A,v0 [,options])
    %
    % options(struct) are
    %   max_iter: maximum number of SCF iterations
    %   tol: tolerance, required maximum norm of residual
    %   k: number of vectors used in DIIS step
    % returns
    %   v: the eigenvector of the last SCF iterations
    %   l: the corresponding eigenvalue
    %   hist: struct containing the residual history
    
    warning('off');
    max_iter = 100;
    tol = 1e-14;
    k = 2;
    if nargin == 3
        if isfield(options, 'max_iter')
            max_iter = options.max_iter;
        end
        if isfield(options, 'tol')
            tol = options.tol;
        end
        if isfield(options, 'k')
            k = options.k;
        end
    end
    
    [v,l] = eigs(A(v0),1,'smallestreal');
    v = v/norm(v);
    res = norm(A(v)*v-l*v);
    ress = A(v)*v-l*v;
    differ = v-v0;
    Vt = v; 
    i = 1;
    while( norm(res(end))>tol && length(res)<max_iter )
        if i<=2
            [v_t,l] = eigs(A(v),1,'smallestreal');
        else
            [v_t,l] = eigs(A(v),1,l);
        end
%         [V,D] = eig(full(A(v)));
%         v_t = V(:,1);
%         l = D(1,1);
        v_t = v_t/v_t(1);
        v_t = v_t/norm(v_t);
        Vt = [Vt v_t];
        differ = [differ v_t-v];
        ress = [ress A(v_t)*v_t-l*v_t];
        if i>=k
            rr = ress(:,end-k+1:end);
            T = [rr'*rr ones(k,1); ones(1,k) 0];
            rhs = [zeros(k,1);1];
            rslt = T\rhs;
            c = rslt(1:end-1);
            v = Vt(:,end-k+1:end)*c;
        else
            v = v_t;
        end
        res = [res norm(A(v)*v-l*v)];
        i = i + 1;
    end
    hist = struct('res',res);
    warning('on');
end