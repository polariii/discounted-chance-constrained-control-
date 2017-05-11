clc;clear all
%% Offline parameters

A=0.5; % whether it is larger than 1 affects the problem
B=1;C=0.5; % system dynamics
nx=size(A,1); % number of states
nu=size(B,2); % number of inputs
N=15;% length of horizon 
W=0.2;% variance of the noise
Gamma=0.9;% discount rate
E0=0.25;% base bound on the output 
t=1; % constraint on the output, P(|y(i)|>=t)
coef=[];
for i=0:N
    coef=[coef,Gamma^i];
end

% important matrices
Q=C'*C;R=1;% weighting matrix
I=[eye(nx), zeros(nx,N*nu)];
K=-dlqr(A,B,Q,R); % nu*nx
E=[eye(nu), zeros(nu,(N-1)*nu)];
O1=zeros(nu,nx);
Qp=I'*Q*I; Rp=[O1 E]'*R*[O1 E];

f=A+B*K; 
M1=zeros(nu*(N-1),nu);M2=eye(nu*(N-1));M3=zeros(nu,nu*N);
M=[M1,M2;M3];
O2=zeros(N*nu,nx);
F=[A  , B*E;
   O2 , M  ];
P=dlyap(F',Qp+Rp);


% % some parameters
r=0.2; % reference value of the output
xbar_r=r/C; % reference value of the mean of the state
q_r=((1-A)*xbar_r)/B; % reference value of the mean of the input
% Kbar=-dlqr(A,B,Q,R);Kn=-dlqr(A,B,Q,R); % feedback gain of second mode
% % solving discrete time Lyapunov equation
% Pbar=dlyap((A+B*Kbar)',Q+Kbar'*R*Kbar); % (A+B*Kbar)'*Pbar*(A+B*Kbar)-Pbar+Q+Kbar'*R*Kb##`    1%% Optimisation: Calculate m*=(m(0),m(1),...,m(N-1))

% initial conditions
x0=-0.8;% initial state
xbar_0=-0.8;% initial mean of the state
X0=0; % 1*eye(Nx,Nx);% initial covariance matrix of the state
Epsilon(1,:)=E0/(1-Gamma);% bound on the output for a inifinite horizon, which will be updated 
% initialising the value of X(i),i=1,...,N-1

X=zeros(N,1);% initialisation of the covariance matrix 
X(1)=(A+B*K)*X0*(A+B*K)'+W;
for i=1:N-1
X(i+1)=(A+B*K)*X(i)*(A+B*K)'+W;
end
for k=1:1
% rng(k) seed for generate random sequence
% Computation
T=50; % time which realisation is up to
for j=1:T
    
cvx_begin quiet 
% cvx_solver sedumi
cvx_precision best
variable bet0 nonnegative 
variable bet(N) nonnegative
variable z(nx+N*nu) % [xbar0-xbar_r;q0-q_r;q1-q_r;q2-q_r;...;qN-1-q_r];
dual variable d1{N}
dual variable d2

ref=[xbar_r;q_r.*ones(N,1)]; % this might be a problem if dimension grows
expression zp(nx+N*nu)
zp=z+ref; % [xbar0;q0;q1;q2;q3;...;qN-1]

expression xbar(N)
for i=1:N
xbar(i)=[A^i,H_matrix(i,A,B)]*zp(1:i+1);
end

minimise ( z'*P*z )  % eigenvalue of P should be nonnegative 

subject to
% 1)
zp(1)==xbar_0;

% 2)
C*X0*C'+(C*xbar_0)^2<=bet0*t^2;
for i=1:N
   d1{i}: C*X(i)*C'+(C*xbar(i))^2<=bet(i)*t^2;
end

% 3)
d2: coef*[bet0;bet]<=Epsilon(j,:)-(Gamma^(N+1))/(1-Gamma);

cvx_end
  
cost(j,:)=cvx_optval;

sign=['The problem is ',cvx_status];
disp(sign)

if string(cvx_status) ==string('Solved')
else
    if string(cvx_status) ==string('Inaccurate/Solved')
    else
    disp('Something went wrong !!!')
    break
    end
end



% update initial conditions
 q0=zp(2);
 w=normrnd(0,sqrt(W)); % noise realisation 
 x0=A*x0+B*K*(x0-xbar_0)+B*q0+w;
 xbar_0=x0;
 
% added term to adjust the mean of control input
 V_noise=K*w;
 for i=1:N-2
     V_noise=[V_noise;K*(f^i)*w]; 
 end
 
 z_new=[xbar_0;zp(nx+nu+1:nx+N*nu)+V_noise;q_r];

 for i=1:N
 xbar_new(i)=[A^i,H_matrix(i,A,B)]*z_new(1:i+1);
end
 
 % update the Epsilon
 e0=(C*X0*C'+(C*xbar_0)^2)/(t^2);
 for i=1:N
 e(i,:)=(C*X(i)*C'+(C*xbar_new(i))^2)/(t^2);
 end
 
 Epsilon(j+1,:)=coef*[e0;e]+(Gamma^(N+1))/(1-Gamma);
 
x(j,:)=x0;

end

figure
plot(Epsilon(1:T),'k')

end

figure
plot(x)

% linsysolve: solution contains NaN or inf