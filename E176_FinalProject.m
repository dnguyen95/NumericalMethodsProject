%% E176 Final Project
% Perturbation Method for Dynamical Systems
% Daniel Nguyen and Austin Chun

% System constants
m1 = 1; m2 = 1;                     % kg
kc1 = 5.8; kc2 = 5.8; kc3 = 5.8;     % N/m
k1 = 5; k2 = 5;                     % N/m
c1 = 0.2; c2 = 0.2;                 % Ns/m

F0 = 1; % Arbitray constant
f1 = F0; f2 = F0;

% GE Matrix form

A = [0,     1,      0,      0;
    -(k1+kc1+kc2)/m1, -c1/m1,   kc2/m1,     0;
    0,      0,      0,      1;
    kc2/m2,    0,   -(k2+kc2+kc3)/m2,      -c2/m2];

B = eye(4);
f = [0; f1; 0; f2];


% Solve for Eigegnvalues/vetors
[U,D] = eig(A);

[~,perm]=sort(diag(D));
D = D(perm,perm);
U = U(:,perm);

[V,~] = eig(A.');
V = V(:,perm);

lam = diag(D);

% Normalize vectors
V = eye(4)/U;

% Forcing function 
Q = V.' * B*f

t = 0:0.01:50;

eta1 = Q(1)/D(1,1) * (1 - exp(lam(1)*t));
eta2 = Q(2)/D(2,2) * (1 - exp(lam(2)*t));
eta3 = Q(3)/D(3,3) * (1 - exp(lam(3)*t));
eta4 = Q(4)/D(4,4) * (1 - exp(lam(4)*t));
eta = [eta1; eta2; eta3; eta4];

x = U*eta;

%% Plotting
figure(1)
plot(t,x(1), t,x(2), t,x(3), t,x(4))
grid on
legend('x1','v1','x2','v2')


%eta = Q./lam .* (1 - exp(



    
