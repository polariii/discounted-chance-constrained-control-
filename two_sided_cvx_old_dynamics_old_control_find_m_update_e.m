clc;clear all
%% Offline parameters

A=2;B=1;C=0.5; % system dynamics
N=15;% length of horizon 
W=0.2;% variance of the noise
Q=C'*C;R=1;% weighting matrix
Gamma=0.9;% discount rate
E0=0.25;% base bound on the output 
t=1; % constraint on the output, P(|y(i)|>=t)
coef=[];
for i=0:N
    coef=[coef,Gamma^i];
end

% some parameters
r=0.0; % reference value of the output
xbar_r=r/C; % reference value of the mean of the state
m_r=((1-A)*xbar_r)/B; % reference value of the mean of the input
Kbar=-dlqr(A,B,Q,R);Kn=-dlqr(A,B,Q,R); % feedback gain of second mode
% solving discrete time Lyapunov equation
Pbar=dlyap((A+B*Kbar)',Q+Kbar'*R*Kbar); % (A+B*Kbar)'*Pbar*(A+B*Kbar)-Pbar+Q+Kbar'*R*Kbar=0
P=dlyap((A+B*Kn)',Q+Kn'*R*Kn); % (A+B*Kn)'*P*(A+B*Kn)-P+Q+Kn'*R*Kn=0

%% Optimisation: Calculate m*=(m(0),m(1),...,m(N-1))

% initial conditions
x0=-0.8;% initial state
xbar_0=-0.8;% initial mean of the state
X0=0; % 1*eye(Nx,Nx);% initial covariance matrix of the state
Epsilon(1,:)=E0/(1-Gamma);% bound on the output for a inifinite horizon, which will be updated 
% initialising the value of X(i),i=1,...,N-1
K=-dlqr(A,B,Q,R);% initial guess for K(i)
f=A+B*K;
X=zeros(N,1);% initialisation of the covariance matrix 
X(1)=(A+B*K)*X0*(A+B*K)'+W;
for i=1:N-1
X(i+1)=(A+B*K)*X(i)*(A+B*K)'+W;
end

% Computation
T=200; % time which realisation is up to
% rng(1) seed for generating random sequence
for j=1:T
    
cvx_begin quiet 
cvx_precision best
variable bet0 nonnegative 
variable bet(N) nonnegative
variable m0
variable m(N-1)
dual variable d1{N}
dual variable d2

expression xbar_n 
expression h(N)
xbar_n=(A^N)*xbar_0+H_matrix(N,A,B)*[m0;m];
h(1)=C*(A*xbar_0+B*m0);
for i=2:1:N
h(i)=C*( (A^i)*xbar_0+H_matrix(i,A,B)*[m0;m(1:i-1)] ) ; 
end

minimise( sum(( h(1:N-1)-r ).^2) + R*(m0-m_r)^2 + sum(R.*( m-m_r ).^2) + (xbar_n-xbar_r)'*Pbar*(xbar_n-xbar_r) )
subject to

% 1)
       C*X0*C'+(C*xbar_0)^2<=bet0*t^2;

for i=1:N
d1{i}: C*X(i)*C'+h(i)^2<=bet(i)*t^2; 
end

% 2)

d2: coef*[bet0;bet]<=Epsilon(j)-(Gamma^(N+1))/(1-Gamma);

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
 w=normrnd(0,sqrt(W));
 x0=A*x0+B*K*(x0-xbar_0)+B*m0+w;
 xbar_0=x0;
 
 % added term to adjust the mean of control input
 V_noise=[];
 for i=1:N-2
     V_noise=[V_noise;K*(f^i)*w]; 
 end
  
 mp0=m(1)+K*w;
 mp=[m(2:N-1)+V_noise;m_r];
 
hp(1)=C*(A*xbar_0+B*mp0);
for i=2:1:N
hp(i)=C*( (A^i)*xbar_0+H_matrix(i,A,B)*[mp0;mp(1:i-1)] ) ; 
end
 
 % update the Epsilon
 e0=(C*X0*C'+(C*xbar_0)^2)/(t^2);
 for i=1:N
 e(i,:)=(C*X(i)*C'+hp(i)^2)/(t^2);
 end
 
 Epsilon(j+1,:)=coef*[e0;e]+(Gamma^(N+1))/(1-Gamma);
 
x(j,:)=x0;

end

%% Analyse the trend of Epsilon

figure
plot(Epsilon(1:T),'r')
title('Epsilon')

figure
plot(cost)
title('Cost')

figure 
plot(x)
title('State Realisation')
