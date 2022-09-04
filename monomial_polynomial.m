function [C0,C1] = monomial_polynomial(AA,ll,vv)
    % monomial_polynomial Build the linearization of a monomial interpolation
    % polynomial for AA(l,v)
    %
    % [C0,C1] = monomial_polynomial(AA,ll,vv)
    %
    % AA(l,v) is the NEPv that needs to be interpolated.
    % ll is a list of approximated eigenvalues
    % vv is a list of corresponding approximated eigenvectors
    A1 = AA(ll(1),vv(:,1));
    A2 = AA(ll(2),vv(:,2));
    C1 = (A2-A1)/(ll(2)-ll(1));
    C0 = C1*ll(1)-A1;
end
