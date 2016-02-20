% TEST - Bulirsch-Stoer Integration method
%
% For method details:
% >> help BulirschStoer
%
% OUTLINE:
%   1) Set up an initial value problem
%   2) Solve using Bulirsch-Stoer
%   3) Solve using ode45
%   4) Compare the solution
%

tSpan = [0,10];   % Time span for initial value problem

z0 = [0.6; 2.5];   %Initial state

% Dynamical system - driven damped pendulum
dynFun = @(t,z)( [z(2,:);  cos(t) - 0.1*z(2,:) - sin(z(1,:))] );

% Solve using Bulirsch-Stoer method
tol = 1e-12;
t = linspace(tSpan(1), tSpan(2), 25);
[z, info] = BulirschStoer(dynFun,t,z0,tol);

% Solve using ode45:
options = odeset('AbsTol',1e-12, 'RelTol',1e-12);
sol = ode45(dynFun,tSpan,z0,options);
tSoln = linspace(tSpan(1), tSpan(2), 150);
zSoln = deval(sol,tSoln);
zCheck = deval(sol,t);

%%%% Plot!

figure(1); clf;

subplot(2,2,1); hold on;
plot(tSoln,zSoln(1,:),'k-')
plot(t,z(1,:),'ro','MarkerSize',10,'LineWidth',2);
legend('ode45','BulirschStoer')
ylabel('angle')
title('Solution')

subplot(2,2,3); hold on;
plot(tSoln,zSoln(2,:),'k-')
plot(t,z(2,:),'ro','MarkerSize',10,'LineWidth',2);
legend('ode45','BulirschStoer')
xlabel('time')
ylabel('rate')


subplot(2,2,2); hold on;
plot(t,abs(z(1,:) - zCheck(1,:)),'ko','MarkerSize',10,'LineWidth',2);
ylabel('angle')
title('| Ode45 - BulirschStoer |')
set(gca,'yScale','log')

subplot(2,2,4); hold on;
plot(t,abs(z(2,:) - zCheck(2,:)),'ko','MarkerSize',10,'LineWidth',2);
xlabel('time')
ylabel('rate')
set(gca,'yScale','log')