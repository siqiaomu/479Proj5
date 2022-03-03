X = usol.'; X1 = X(:,1:end-1); X2 = X(:,2:end);


[U, Sigma, V] = svd(X1,'econ');
S = U'*X2*V*diag(1./diag(Sigma));
[eV,D] = eig(S);
mu = diag(D);
omega = log(mu)/(dt);
Phi = U*eV;


y0 = Phi\u; %psuedo-inverse initial conditions

u_modes = zeros(size(y0,1),length(t));
for iter = 1:length(t)
    u_modes(:,iter) = (y0.*exp(omega*t(iter)));
end


u_dmd = Phi*u_modes;