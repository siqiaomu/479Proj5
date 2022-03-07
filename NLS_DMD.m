% L is the spread of the x points we are looking at (-20,20)
% n gives the discrete number of points along that interval we look at
% slices gives number of time snapshots (after zero) we will simulate 

% note, for whatever reason, that this is the number of time snapshots on 
% the interval of (0, 2*pi) in t
L = 40; n = 512; slices = 20;

time_int = [0 2*pi];

% generate data (usol) with initial condition u, at time points t,
% with separation dt
% note: the xx and tt outputs are just helpful for plotting the resulting
% plot with the surf function
[t,usol,u,dt,tt,xx] = nls_data(L,n,slices,time_int);

% set up the matrices for DMD
X = usol.'; X1 = X(:,1:end-1); X2 = X(:,2:end);

% perform singular value decomposition on X1
[U, Sigma, V] = svd(X1,'econ');

% get tilde S
S = U'*X2*V*diag(1./diag(Sigma));

% get eigenvalues of tilde S
[eV,D] = eig(S); % eV contains the 'right eigenvectors' 
% and D a diagonal matrix of the generalized eigenvalues
mu = diag(D); % get the eigenvalues out of D

omega = log(mu)/(dt); % 'transform' eigenvalues before constructing dmd modes
Phi = U*eV; % calculate DMD modes

y0 = Phi\u; % use psuedo-inverse to obtain initial conditions

% construct the predicted dynamics (at the same time snapshots we
% simulated our data for)
u_modes = zeros(size(y0,1),length(t));
for iter = 1:length(t)
    u_modes(:,iter) = (y0.*exp(omega*t(iter)));
end

u_dmd = Phi*u_modes;

% some of the following is optional (related to the plots I made)

% determine (indices of) which eigenvalues are growth modes/ which decay
%abs_mu = abs(mu);
%grow = find(abs_mu >= 1);
%decay = find(abs_mu < 1);

%nls_plt_modes;

% let's look at the error of the solution, especially when we predict out
% further in time

slices = 350; time_int = [0 100];

% exact simulation
[t_2,usol_2,u_2,dt_2,tt_2,xx_2] = nls_data(L,n,slices,time_int);

% get new ICs
y0_2 = Phi\u_2;


% now use the previously computed dmd approximation to simulate on this
% larger time interval

u_modes_2 = zeros(size(y0,1),length(t_2));
for iter = 1:length(t_2)
    u_modes_2(:,iter) = (y0.*exp(omega*t_2(iter)));
end

u_dmd_2 = Phi*u_modes_2;

% compute the norm difference
dif = u_dmd_2' - usol_2;
for iter = 1:length(t_2)
    norm_dif(iter) = norm(dif(iter,:));
end

abs_dif = abs(dif);

figure(1)
surf(tt_2,xx_2,log10(abs_dif'));

figure(2)
plot(tt_2,norm_dif,'linewidth',2.5,'Color','black');



