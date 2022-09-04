function [C0,C1] = newton_polynomial(AA,ll,vv,throwError)
    % newton_polynomial Build the linearization of a Newton interpolation
    % polynomial for AA(l,v)
    %
    % [C0,C1] = newton_polynomial(AA,ll,vv)
    %
    % AA(l,v) is the NEPv that needs to be interpolated.
    % ll is a list of approximated eigenvalues
    % vv is a list of corresponding approximated eigenvectors
    
    [n,N] = size(vv);
    C0 = cell(N-1,N-1);
    C0(:) = {zeros(n)};
    C1 = cell(N-1,N-1);
    C1(:) = {zeros(n)};
    if nargin==3
        throwError = false;
    end
    Ai = divided_differences(AA,ll,vv,throwError);
    for i = 1:N-1
        C0{1,i} = Ai{i,1};
        if i<N-1
            C0{i+1,i} = ll(i)*eye(n);
            C0{i+1,i+1} = eye(n);
            C1{i+1,i} = eye(n);
        end
    end
    C0{1,end} = C0{1,end} - ll(end-1)*Ai{end,1};
    C1{1,end} = -Ai{end,1};
    C0 = cell2mat(C0);
    C1 = cell2mat(C1);
end

function [Ai] = divided_differences(AA,ll,vv,throwError)
    [n,N] = size(vv);
    Y = zeros(n*N,n);
    for i = 1:length(ll)
        Y((i-1)*n+1:i*n,:) = AA(ll(i),vv(:,i));
    end
    M = zeros(N);
    M(:,1) = ones(N,1);
    for j = 2:N
        M(j:end,j) = ll(j:end)-ll(j-1);
    end
    M = cumprod(M,2);
    if throwError && rcond(M)<1e-15
       ME = MException(1,'Divided differences has a singular coefficient matrix');
       throw(ME);
    end
    
    M = kron(M,eye(n));
    Ai = M\Y;
    Ai = mat2cell(Ai,n*ones(1,N),n);
end