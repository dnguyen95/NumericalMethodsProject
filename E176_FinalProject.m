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

% Define number of STATES in the system
N = 4;

% GE Matrix form
A = [0,     1,      0,      0;
    -(k1+kc1+kc2)/m1, -c1/m1,   kc2/m1,     0;
    0,      0,      0,      1;
    kc2/m2,    0,   -(k2+kc2+kc3)/m2,      -c2/m2];
B = eye(4);
f = [0; f1; 0; f2];


% Solve for Eigegnvalues/vetors
[U,D] = eig(A);
% Sort eigenvalues/vectors for 
[~,perm]=sort(diag(D));
D = D(perm,perm);
U = U(:,perm);

[V,~] = eig(A.');
V = V(:,perm);

lam = diag(D);

% Normalize vectors
Vt = eye(4)/U;
V = Vt.';

% Forcing function 
Q = V.' * B*f;

t = 0:0.1:50;
eta = zeros(N, length(t));
for i = 1:N
   eta(i,:) = Q(i)/lam(i) * (1 - exp(lam(i)*t)); 
end

%{
eta1 = Q(1)/lam(1) * (1 - exp(lam(1)*t));
eta2 = Q(2)/lam(2) * (1 - exp(lam(2)*t));
eta3 = Q(3)/lam(3) * (1 - exp(lam(3)*t));
eta4 = Q(4)/lam(4) * (1 - exp(lam(4)*t));
eta = [eta1; eta2; eta3; eta4];
%}

%% Solve for x from decouple eta's
x = U*eta;

%% Plotting
figure(1)
plot(t,x(1,:), t,x(2,:), t,x(3,:),'--', t,x(4,:),'--')
grid on
legend('x_1(t)','v_1(t)','x_2(t)','v_2(t)')


%%%%%%%%%%%%%%%%%%%
%% Perturbation %%
%%%%%%%%%%%%%%%%%%%
% Deviations in system parameters
dm1 = 0.3229; dm2 = 0.2253; dkc1 = -0.1556;
dkc2 = 0.0917; dkc3 = -0.1843; dk1 = -0.1448;
dk2 = 0.0531; dc1 = 0.0199; dc2 = -0.1476;

% Recalulate system parameters
m1=m1+dm1; m2=m2+dm2; kc1=kc1+dkc1;
kc2=kc2+dkc2; kc3=kc3+dkc3; k1=k1+dk1;
k2=k2+dk2; c1=c1+dc1; c2=c2+dc2;

% Calculate new A
A_new = [0,     1,      0,      0;
    -(k1+kc1+kc2)/m1, -c1/m1,   kc2/m1,     0;
    0,      0,      0,      1;
    kc2/m2,    0,   -(k2+kc2+kc3)/m2,      -c2/m2];
% Determine the deviation matrix
dA = A_new - A;

%% Calculate deviations in eigenvalues/vectors




    
