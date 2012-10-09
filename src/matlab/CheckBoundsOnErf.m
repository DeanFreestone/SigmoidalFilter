%% here I will try to make a bound on one of the error function at the
%% start of the derivation.

clear
close all
clc

x = -20:0.001:20;

Q = 0.5 - 0.5*erf(x/sqrt(2));   % for x from 0 to infty 

g = 1 - Q;          % this is our g(x) for x=-infty to 0 g(x) <= 0.5*exp(-x.^2/2)

bound = 0.5*exp(-x.^2/2);

bound2 = 1 - bound;

% l_bound = x.*exp(-x.^2/2)./((1+x.^2)*sqrt(2*pi));
l_bound = 1.*exp(-x.^2/2)./((1+x.^2)*sqrt(2*pi));

l_bound = 1- l_bound;

figure
% plot(x,Q,x,bound,x,g,x,bound2,x,l_bound)
% legend('Q','bound','g','bound 2','lower bound')
plot(x,bound,x,g,x,l_bound)
legend('bound','g','lower bound')


%% 
beta = 1.08;
l_bound2 = exp(-beta*x.^2)*sqrt(2*exp(1)/pi)*sqrt(beta-1)/beta;
l_bound2 = 1-l_bound2;

figure
% plot(x,Q,x,bound,x,g,x,bound2,x,l_bound)
% legend('Q','bound','g','bound 2','lower bound')
plot(x,bound,x,g,x,l_bound2)
legend('bound','g','lower bound chernoff')