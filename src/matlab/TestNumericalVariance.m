
clc
clear
close all

min_w = 0;
dw = 0.001;
max_w = 20;

min_z = 0;
dz = 0.001;
max_z = 20;

w = min_w:dw:max_w;
z = min_z:dz:max_z;

mu = 1;
s = 2;

%% this is a the point before bounding the variance
coeff = exp(-mu^2/(2*s^2+1))/(2*pi*sqrt(2*s^2+1));

just_z_bit = exp((2*z*mu - (s^2 + 1)*z.^2)/(4*s^2+2));

just_w_bit = exp((2*mu*w - (s^2+1)*w.^2)/ (4*s^2+2));

zw_coeff = 2*s^2*w / (4*s^2+2);

plot(z,just_z_bit)

temp = zeros(1,length(z));          % pre-allocate for speed
for n=1:length(z)
    temp(n) = sum( just_w_bit.*exp(zw_coeff*z(n)) ) * dw;
end

E_y2_int_1 = sum(coeff*just_z_bit.*temp) * dz;

%% this bit has epsilon in it
epsilon = 1;
coeff = exp(-mu^2/(2*s^2+1))/(2*pi*sqrt(2*s^2+1));

just_z_bit = exp((2*mu*z - (s^2 + 1)*z.^2 + s^2*z.^2/epsilon^2) / (4*s^2+2));

just_w_bit = exp((2*mu*w - (s^2+1)*w.^2 + s^2*w.^2*epsilon^2) / (4*s^2+2));

int_z = sum(just_z_bit)*dz;
int_w = sum(just_w_bit)*dw;
E_y2_int_2 = coeff*int_z*int_w;

%% this at the first step after integrating out x
coeff = exp(-(mu^2)/(2*s^2))/(2*pi*sqrt(2*s^2+1));
temp = zeros(1,length(z));          % pre-allocate for speed

for n=1:length(z)
    
    temp(n) = sum(exp(- 0.5*(w.^2 + z(n)^2) + (s^2*w + s^2*z(n) + mu).^2 / (s^2*(4*s^2+2))));

end

E_y2_int = coeff*sum(temp)*dw*dz;

%% afer integrating out w and x
% 
% coeff  = exp(-mu^2/(2*s^2+1))/(2*pi*sqrt(2*s^2+1));
% 
% firstbit = exp((2*mu*z - (s^2+1)*z.^2)/(4*s^2+2) + (mu + s^2*z).^2/(2*(s^2+1)*(2*s^2+1)));
% secondbit = sqrt(pi*(1-1/(2*s^2+2))) - sqrt(pi/2)*sqrt(2*s^4+3*s^2+1)/(s^2+1)*erf((-z*s^2-mu)/sqrt(4*s^4+6*s^2+2));
% 
% figure,plot(firstbit.*secondbit)
% 
% E_y_simp = coeff*sum(firstbit.*secondbit)*dz;

% %% integrated out x and w and bounding erf
% 
% erf_argument = z;%(-z*s^2-mu)/sqrt(4*s^4+6*s^2+2);
% temp = 1-erf(erf_argument);
% temp2 = erfc(erf_argument);
% temp3 = exp(-erf_argument.^2);
% 
% figure,plot(z,temp,z,temp2,z,temp3)
% legend('from erf','from erfc','from bound')

%% using the approximation of the error function
dx = 0.001;
x_neg = -20:dx:0;
x_pos = 0:dx:20;

x = -20:dx:20;
g = 0.5+0.5*erf(x/sqrt(2));
f = exp(-(x-mu).^2/(2*s^2)) / (sqrt(2*pi)*s);

b_u = 1.08;
a_u = sqrt(2*exp(1)/pi)*sqrt(b_u-1)/b_u;

a_l = 0.5;
b_l = 0.5;

f_pos = exp(-(x_pos-mu).^2/(2*s^2)) / (sqrt(2*pi)*s);
f_neg = exp(-(x_neg-mu).^2/(2*s^2)) / (sqrt(2*pi)*s);

g_u = 1-a_u*exp(-b_u*(x_pos).^2);
g_l = a_l*exp(-b_l*(x_neg).^2);

figure
plot(x_neg,f_neg,x_neg,g_l, x_pos,f_pos,x_pos,g_u,x,g)
legend('pdf from -inf to 0','erf from -inf to 0', 'pdf from 0 to inf', 'erf from 0 to inf','sigmoid')

temp1 = g_l.^2.*f_neg;
temp2 = g_u.^2.*f_pos;
% temp3 = g.^2.*f;
% figure
% plot(x_neg,temp1,x_pos,temp2,x,temp3)
% 
% sum(temp3)*dx
E_y_erf_bound = sum(temp1)*dx+sum(temp2)*dx;

%%

N_realisations = 1e8;
x = mu + s*randn(1,N_realisations);
y = 0.5*(erf(x/sqrt(2)) + 1);           % transform RV through error function
mean(y.^2)
numerical_M = mean(y);
numerical_V = var(y);
disp(['Numerical Mean (Monte Carlo) = ' num2str(numerical_M)])
disp(['Numerical Variance (Monte Carlo) = ' num2str(numerical_V)])

%% Now test the Unscented transform
% create sigma points
n=1;        % number of states
X = [mu + sqrt(n)*s, mu - sqrt(n)*s];

% pass sigma points through nonlinearity
g_X = 0.5*(erf(X/sqrt(2)) + 1);

% approximate the mean
E_X = (1/(2*n))*sum(g_X);
disp(['The expected value from simple UT is ' num2str(E_X)])

% approximate the covariance
P_X = sum((g_X - E_X)*(g_X - E_X)') / (2*n);
disp(['The covariance from simple UT is ' num2str(P_X)])

%%

z = mu/sqrt(1+s^2);
M = 0.5*(erf(z/sqrt(2)) + 1);
disp(['Analytic Mean (Frey) = ' num2str(M)])

V_temp = mu/sqrt(1+s^2);
Phi_V_temp = 0.5*(erf(V_temp/sqrt(2)) + 1);
V = Phi_V_temp*(1-Phi_V_temp)*s^2/(s^2 + pi/2);
disp(['Analytic Upper Bound of Variance (Frey) = ' num2str(V)])

% E_y2 = erfc(mu/sqrt(4*s^2+2))^2*(2*s^2 + 1)/(4*sqrt(2*s^2+1));
% % 
% V2 = E_y2 - M^2;
% disp(['My Analytic Variance  = ' num2str(V2)])

V3 = E_y2_int - M^2;
disp(['Numerically Integrated Variance  = ' num2str(V3)])

V4 = E_y2_int_1 - M^2;
disp(['Numerically Integrated Variance - after simplification = ' num2str(V4)])
%%
V5 = E_y2_int_2 - M^2;
disp(['Numerically Integrated Variance - after epsilon = ' num2str(V5)])

%%
% V6 = E_y_simp - M^2;
% disp(['Numerically Integrated Variance - reduced to z = ' num2str(V6)])

%%
V7 = E_y_erf_bound - M^2;
disp(['Numerically Integrated bound on erf = ' num2str(V7)])

%%
% clc
% clear
% close all
% 
% % test bound on error function
% x = linspace(0,10,100000);
% % error_func = erf(x/sqrt(2));
% % upr_bnd = -2*exp(-x.^2/2)./(x*sqrt(2*pi)) + 1;
% 
% Q = 0.5 - 0.5*erf(x/sqrt(2));
% bound_Q = exp(-x.^2/2)./(x*sqrt(2*pi));
% display(['Q function sum: ' num2str(sum(Q(2:end)))])
% display(['bound Q function: ' num2str(sum(bound_Q(2:end)))])
% 
% Chernoff_bound = 0.5*exp(-x.^2/2);
% display(['Chernoff Q sum: ' num2str(sum(Chernoff_bound(2:end)))])
% 
% figure
% plot(x,Q,x,bound_Q,x,Chernoff_bound)
% legend('Q function','Q bound','Chernoff')

