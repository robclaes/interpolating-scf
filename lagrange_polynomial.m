function [C0,C1] = lagrange_polynomial(AA,ll,vv)
    % newton_polynomial Build the linearization of a Lagrange interpolation
    % polynomial for AA(l,v)
    %
    % [C0,C1] = lagrange_polynomial(AA,ll,vv)
    %
    % AA(l,v) is the NEPv that needs to be interpolated.
    % ll is a list of approximated eigenvalues
    % vv is a list of corresponding approximated eigenvectors
    
    [n,N] = size(vv);
    C0 = cell(N-1,N-1);
    C0(:) = {zeros(n)};
    C1 = cell(N-1,N-1);
    C1(:) = {zeros(n)};
    Ai = @(i) AA(ll(i),vv(:,i));
    for i = 1:N-1
        C0{1,i} = ll(i+1)*Ai(i);
        C1{1,i} = Ai(i);
        if i<N-1
            C0{i+1,i} = ll(i)*eye(n);
            C0{i+1,i+1} = -ll(i+2)*theta(i+1,ll)*eye(n);
            C1{i+1,i} = eye(n);
            C1{i+1,i+1} = -theta(i+1,ll)*eye(n);
        end
    end
    C0{1,end} = C0{1,end} + ll(end-1)*Ai(N)/theta(N,ll);
    C1{1,end} = C1{1,end}+Ai(N)/theta(N,ll);
    C0 = cell2mat(C0);
    C1 = cell2mat(C1);
end

function t = theta(i,ll)
    lli = ll(i)-ll;
    lli_m1 = ll(i-1)-ll;
    lli(i) = 1;
    lli_m1(i-1) = 1;
    wi = prod(lli);
    wi_m1 = prod(lli_m1);
    t =wi/wi_m1;
end