function [Af] = gen_gpe(f, N, l, omega)
    xn = linspace(-l,l,N+2);
    xn = xn(2:end-1);
    yn = xn;
    h = (2*l)/(N+1);
    e = ones(N,1);
    DN = spdiags([-1/2*e 0*e 1/2*e],-1:1,N,N);
    D2N = spdiags([e -2*e e],-1:1,N,N);
    
    M = kron(D2N, eye(N))+kron(eye(N),D2N);
    Ms = h*kron(diag(yn),DN)-h*kron(DN, diag(xn));
    [X,Y] = meshgrid(xn,yn);
    F = f(X,Y).';
    ff = h^2*F(:);
    Af = spdiags(ff,0,length(ff),length(ff)) - 1/2*M - 1i*omega*Ms;
end