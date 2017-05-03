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

[V,D] = eig(A.');
%[~,perm]=sort(diag(D));
V = V(:,perm);

lam = diag(D);

% Normalize vectors
VtU = V.'*U;
for i = 1:N
    U(:,i) = U(:,i) / sqrt(VtU(i,i));
    V(:,i) = V(:,i) / sqrt(VtU(i,i));
end

% Forcing function 
Q = V.' * B*f;

t = 0:0.1:50;
eta = zeros(N, length(t));
for i = 1:N
   eta(i,:) = Q(i)/lam(i) * (1 - exp(lam(i)*t)); 
end

% Solve for x from decouple eta's
x = U*eta;

% Plotting
figure(1)
%subplot(3,1,1)
plot(t,x(1,:), t,x(2,:), t,x(3,:),'--', t,x(4,:),'--')
xlabel('Time [ s ]')
ylabel('Amplitude [ m or m/s ]')
title('System Forced Response Exact, no perturbation')
grid on
legend('x_1(t)','v_1(t)','x_2(t)','v_2(t)')
set(gcf,'color','white')


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
A_o = A; % Label old one as old
A = [0,     1,      0,      0;
    -(k1+kc1+kc2)/m1, -c1/m1,   kc2/m1,     0;
    0,      0,      0,      1;
    kc2/m2,    0,   -(k2+kc2+kc3)/m2,      -c2/m2];
% Determine the deviation matrix
dA = A - A_o;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perturbation Analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate deviations in eigenvalues/vectors
lam_o = lam;
U_o = U;
V_o = V;

% New Eigenvalues
lam = zeros(1,N);
for i = 1:N
   lam(i) = lam_o(i) + V_o(:,i).' * dA * U_o(:,i);
end

% New Eigenvectors
U = zeros(N,N);
for i = 1:N
    dU = 0;
    for k = 1:N
        if(k ~= i)
            dU = dU + ( V_o(:,k).' * dA * U_o(:,i)) / (lam_o(i) - lam_o(k)) * U_o(:,k); 
        end
    end
    U(:,i) = U_o(:,i) + dU;
end

V = zeros(N,N);
for i = 1:N
    dV = 0;
    for k = 1:N
        if(k ~= i)
            dV = dV + ( V_o(:,i).' * dA * U_o(:,k)) / (lam_o(i) - lam_o(k)) * V_o(:,k); 
        end
    end
    V(:,i) = V_o(:,i) + dV;
end

% Normalize vectors
VtU = V.'*U;
for i = 1:N
    U(:,i) = U(:,i) / sqrt(VtU(i,i));
    V(:,i) = V(:,i) / sqrt(VtU(i,i));
end

% Extract new, more precise eigenvalues
D = V.' * A * U;
lam = diag(D);

%% Recalculate system response

% Forcing function 
Q = V.' * B*f;

t = 0:0.1:50;
eta = zeros(N, length(t));
for i = 1:N
   eta(i,:) = Q(i)/lam(i) * (1 - exp(lam(i)*t)); 
end

% Solve for x from decouple eta's
x = U*eta;

% Plotting
figure(2)
%subplot(3,1,2)
plot(t,x(1,:), t,x(2,:), t,x(3,:),'--', t,x(4,:),'--')
xlabel('Time [ s ]')
ylabel('Amplitude [ m or m/s ]')
title('System Forced Response with Perturbation Analysis')
grid on
legend('x_1(t)','v_1(t)','x_2(t)','v_2(t)')
set(gcf,'color','white')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exact Solution to Perturbed System %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve for Eigegnvalues/vetors
[U_exact,D_exact] = eig(A);
% Sort eigenvalues/vectors for 
[~,perm]=sort(diag(D_exact));
D_exact = D_exact(perm,perm);
U_exact = U_exact(:,perm);

[V_exact,D] = eig(A.');
[~,perm]=sort(diag(D));
V_exact = V_exact(:,perm);

lam_exact = diag(D_exact);



% Normalize vectors
VtU = V_exact.'*U_exact;
for i = 1:N
    U_exact(:,i) = U_exact(:,i) / sqrt(VtU(i,i));
    V_exact(:,i) = V_exact(:,i) / sqrt(VtU(i,i));
end


% Forcing function 
Q = V_exact.' * B*f;

t = 0:0.1:50;
eta = zeros(N, length(t));
for i = 1:N
   eta(i,:) = Q(i)/lam_exact(i) * (1 - exp(lam_exact(i)*t)); 
end

% Solve for x from decouple eta's
x = U_exact*eta;

% Plotting
figure(3)
%subplot(3,1,3)
plot(t,x(1,:), t,x(2,:), t,x(3,:),'--', t,x(4,:),'--')
xlabel('Time [ s ]')
ylabel('Amplitude [ m or m/s ]')
title('System Forced Response Exact, with perturbation')
grid on
set(gcf,'color','white')







