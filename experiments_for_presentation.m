clear;
close all;
rng(0);
nb_times = 10;
%% Problem setup

N = 25;
n=N*N;
f = @(x,y) (x.^2+y.^2)/2;
Af = gen_gpe(f,N,1,0.85);
beta = 6;
A = @(v) sparse(Af+ beta*diag(abs(v).^2));
[V,~] = eigs(Af,1,'smallestreal');
v0 = orth(V(:,1));   


%% Standard SCF example

options.tol      = 1e-14;
options.max_iter = 20;
tic();
for j = 1:nb_times
    [v, l, hist]     = SCF(A, v0, options);
end
t = toc();
t_scf = t/nb_times/length(hist.res)
figure;
semilogy(hist.res,"k--"); hold on;

%% Secant SCF example

options.tol        = 1e-14;
options.max_iter   = 20;
options.iter_add   = 1;
options.p_hist     = 2;
options.ploteigs   = false;
options.prev_shift = true;
tic();
for j = 1:nb_times
    [v_s,l_s,hist_s]   = SCF_newton(A, v0, options);
end
t = toc();
t_newton = t/nb_times/length(hist_s.res)
semilogy(hist_s.res,"k-"); 

%% Lagrange Secant SCF example

options.tol        = 1e-14;
options.max_iter   = 9;
options.iter_add   = 1;
options.p_hist     = 2;
options.ploteigs   = false;
options.prev_shift = true;
tic();
for j = 1:nb_times
    [v_s,l_s,hist_s]   = SCF_lagrange(A, v0, options);
end
t = toc();
t_lagrange = t/nb_times/length(hist_s.res)
semilogy(hist_s.res,"b-"); 


%% Monomial Secant SCF example

options.tol        = 1e-14;
options.max_iter   = 20;
options.iter_add   = 1;
options.p_hist     = 2;
options.ploteigs   = false;
options.prev_shift = true;
tic();
for j = 1:nb_times
    [v_s,l_s,hist_s]   = SCF_monomial(A, v0, options);
end
t = toc();
t_monomial = t/nb_times/length(hist_s.res)
semilogy(hist_s.res,"r-"); 

%% DIIS SCF example
options.k = 2;
options.max_iter = 18;
tic();
for j = 1:nb_times
    [v_diis,l_diis,hist_diis] = SCF_DIIS(A,v0,options);
end25
t = toc();
t_DIIS = t/nb_times/length(hist_diis.res)
semilogy(hist_diis.res,"k*-");
legend("Standard SCF", "Secant SCF","Secant SCF lagrange","Secant SCF monomial", "SCF + DIIS");
xlabel("Iteration i");
ylabel("Residual");
