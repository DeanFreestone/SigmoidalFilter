clc
clear
close all

FS_Label = 8;
FS_Tick = 8;

%%
% set the parameters for the sigmoid and the Gaussian RV
mu_x = 0;         % mean of the Gaussian
sigma_x = 9;       % standard deviation of the Gaussian
%%

varsigma = 0.56;    % slope parameter for the sigmoid
v_0 = 0;                % threshold parameter for the sigmoid

%% plot the sigmoid and the distributions

% create the sigmoid
xmin = mu_x-40;
xmax = mu_x+40;
x = linspace(xmin,xmax,10000);
g = 1./(1+exp(varsigma*(v_0 - x)));

% create the Gaussian
f_X_x = exp(-(x - mu_x).^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));
% test = gaussmf(x,[sigma_x,mu_x]) / (sigma_x * sqrt(2*pi));
% create the transformed Gaussian - analytically
y = linspace(1e-10,1-1e-10,1e7);
h_y = -(1/varsigma)*log((1./y) - 1) + v_0;
h_dash_y = -1./(varsigma*y.*(y-1));
f_Y_y = 2* h_dash_y.* exp(-(h_y - mu_x).^2 / (2*sigma_x)) / (sigma_x*sqrt(2*pi));

figure('color','white','units','centimeters','position',[2 2 9 9],'papersize',[9 9],'filename','Mapping.pdf')
subplot(223)
plot(x,g,'k')
ylim([-0.2 1.2])
xlabel('$x$','fontsize',FS_Label,'interpreter','latex')
ylabel('$y$','fontsize',FS_Label,'interpreter','latex')
set(gca,'fontsize',FS_Tick)
xlim([xmin xmax])
box off
axis square

subplot(221)
plot(x,f_X_x,'k')

ylabel('$f_X(x)$','fontsize',FS_Label,'interpreter','latex')
set(gca,'fontsize',FS_Tick,'yticklabel',{})
xlim([xmin xmax])
box off
axis square

subplot(224)
plot(f_Y_y,y,'k')
xlabel('$f_Y(y)$','fontsize',FS_Label,'interpreter','latex')
ylim([-0.2 1.2])
set(gca,'fontsize',FS_Tick,'xticklabel',{})
box off
axis square

hold on

%% check Taylor series for a single realisation

x = mu_x + sigma_x*randn(1,1);                  % generate 1 realisation of the GRV
g = 1./(1+exp(varsigma*(v_0 - x)));              % transform the GRV through sigmoid
disp(['The actual value from propogating one realisation through nonlinearity is ' num2str(g)])

% now find the Taylor series approximation (as below)
% tilde g = g(mu_x) + g'(mu_x)(x-mu_x) + 1/2! g''(mu_x)(x-mu_x)^2 + 1/3!g'''(mu_x)(x-mu_x)^3

g_mu = 1./(1+exp(varsigma*(v_0 - mu_x)));
x_tilde = x - mu_x;
g1 = varsigma*(g_mu - g_mu^2);
g2 = varsigma*g1*(1-2*g_mu);
g3 = varsigma^2*g1*(1-6*g_mu+6*g_mu^2);
g4 = varsigma^3*g1*(1-14*g_mu+36*g_mu^2 - 24*g_mu^3);

tilde_g = g_mu + g1*x_tilde + 0.5*g2*x_tilde.^2 + (1/6)*g3*x_tilde.^3 + (1/24)*g4*x_tilde.^4;
disp(['The approximation (4th-order) from propogating one realisation through nonlinearity is ' num2str(tilde_g)])


%%
% In this part of the code I look at how well the Taylor series
% approximates the sigmoid

x = linspace(xmin,xmax,100000);
g = 1./(1+exp(varsigma*(v_0 - x)));

g_mu = 1./(1+exp(varsigma*(v_0 - mu_x)));
x_tilde = x - mu_x;
g1 = varsigma*(g_mu - g_mu^2);
g2 = varsigma*g1*(1-2*g_mu);
g3 = varsigma^2*g1*(1-6*g_mu+6*g_mu^2);
g4 = varsigma^3*g1*(1-14*g_mu+36*g_mu^2 - 24*g_mu^3);

tilde_g = g_mu + g1*x_tilde + 0.5*g2*x_tilde.^2 + (1/6)*g3*x_tilde.^3 + (1/24)*g4*x_tilde.^4;
tilde_g2 = g_mu + g1*x_tilde + 0.5*g2*x_tilde.^2;

% figure
% plot(x,g,x,tilde_g2)
% ylim([-.5 1.5])

%% Monte Carlo
% check expectation
NRealisations = 10000000;
x = mu_x + sigma_x*randn(1,NRealisations);

g = 1./(1+exp(varsigma*(v_0 - x)));

E_g = mean(g);                   % numerical ('true') expectation
P_g = cov(g);
disp(['The expected value from Monte Carlo simulation is ' num2str(E_g)])
disp(['The covariance from Monte Carlo simulation is ' num2str(P_g)])


% ~~~~~~~~~~~~~~
g_mu = 1./(1+exp(varsigma*(v_0 - mu_x)));
x_tilde = x - mu_x;
g1 = varsigma*(g_mu - g_mu^2);
g2 = varsigma*g1*(1-2*g_mu);
g3 = varsigma^2*g1*(1-6*g_mu+6*g_mu^2);
g4 = varsigma^3*g1*(1-14*g_mu+36*g_mu^2 - 24*g_mu^3);

E_tilde_g = g_mu + g1*mean(x_tilde) + 0.5*g2*mean(x_tilde.^2) + (1/6)*g3*mean(x_tilde.^3) + (1/24)*g4*mean(x_tilde.^4);
disp(['The expected value from forth-order Taylor series Monte Carlo simulation is ' num2str(E_tilde_g)])

E_tilde_g2 = g_mu + g1*mean(x_tilde) + 0.5*g2*mean(x_tilde.^2);
disp(['The expected value from second-order Taylor series Monte Carlo simulation is ' num2str(E_tilde_g)])

%% Now test the Unscented transform
% create sigma points
n=1;        % number of states
X = [mu_x + sqrt(n)*sigma_x, mu_x - sqrt(n)*sigma_x];

% pass sigma points through nonlinearity
g_X = 1./(1+exp(varsigma*(v_0 - X)));

% approximate the mean
E_X = (1/(2*n))*sum(g_X);
disp(['The expected value from simple UT is ' num2str(E_X)])

% approximate the covariance
P_X = sum((g_X - E_X)*(g_X - E_X)') / (2*n);
disp(['The covariance from simple UT is ' num2str(P_X)])


%% create a Gaussian Approximation from the UT
% create the Gaussian
y2 = linspace(-1.2,1.2,1e7);
f_Y_y_UT = exp(-(y2 - E_X).^2 / (2*P_X)) / (sqrt(P_X*2*pi));

plot(f_Y_y_UT,y2,'r')

leg = legend('Actual','UT Approx.');
set(leg,'fontsize',FS_Label,'box','off','units','centimeters','position',[6,4,1,1])
%% General Unscented transform
% kappa = 3-n;
% alpha = 10e-3;
% lambda = alpha^2*(n+kappa)-n;       % need to check this guy
% beta = 2;
% 
% % create sigma points
% X = [mu_x mu_x + sqrt(n + lambda)*sigma_x, mu_x - sqrt(n + lambda)*sigma_x];
% 
% % create the weights
% Wm = [lambda/(n+lambda) 1/(2*(n+lambda)) 1/(2*(n+lambda))];
% Wc = [lambda/(n+lambda)+(1-alpha^2+beta), 1/(2*(n+lambda)), 1/(2*(n+lambda))];
% 
% % pass sigma points through nonlinearity
% g_X = 1./(1+exp(varsigma*(v_0 - X)));
% 
% % approximate the mean
% E_WX = sum(Wm.*g_X);
% disp(['The expected value from generalised UT is ' num2str(E_WX)])
% 
% % approximate the covariance
% P_WX = sum(Wc.*(g_X - E_WX)*(g_X - E_WX)');
% disp(['The covariance from generalised UT is ' num2str(P_WX)])

%%