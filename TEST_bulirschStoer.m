% TEST - Bulirsch-Stoer Integration method
%

tSpan = [0,50];   % Time span

z0 = [0.3; 2.5];   %Initial state

% Dynamical system - driven damped pendulum
dynFun = @(t,z)( [z(2,:);  cos(t) - 0.1*z(2,:) - sin(z(1,:))] );

% Test the bulirsch-stoer method
tol = 1e-14;
t = linspace(tSpan(1), tSpan(2), 20);
[z, info] = BulirschStoer(dynFun,t,z0,tol);

% Accurately solve using ode45:
options = odeset('AbsTol',1e-12, 'RelTol',1e-12);
sol = ode45(dynFun,tSpan,z0,options);
zSoln = deval(sol,t);

%%%% Plot!

figure(1); clf;

subplot(2,2,1); hold on;
plot(t,zSoln(1,:),'k-')
plot(t,z(1,:),'ro','MarkerSize',10,'LineWidth',2);
legend('ode45','BulirschStoer')
ylabel('angle')
title('Solution')

subplot(2,2,3); hold on;
plot(t,zSoln(2,:),'k-')
plot(t,z(2,:),'ro','MarkerSize',10,'LineWidth',2);
legend('ode45','BulirschStoer')
xlabel('time')
ylabel('rate')


subplot(2,2,2); hold on;
plot(t,abs(z(1,:) - zSoln(1,:)),'ko','MarkerSize',10,'LineWidth',2);
ylabel('angle')
title('Error')
set(gca,'yScale','log')

subplot(2,2,4); hold on;
plot(t,abs(z(2,:) - zSoln(2,:)),'ko','MarkerSize',10,'LineWidth',2);
xlabel('time')
ylabel('rate')
set(gca,'yScale','log')