clc
clear
close all

varsigma = 7;
v_0 = 3;

dv = 0.001;
v_neg = -20:dv:v_0;
v_pos = v_0:dv:20;

v = -20:dv:20;

z = (v-v_0)/(varsigma*sqrt(2));
g = 0.5*(erf(z)+1);

b2 = 1.08;
a2 = sqrt(2*exp(1)/pi)*sqrt(b2-1)/b2;

a1 = 0.5;
b1 = 0.5;
z_neg = (v_neg-v_0)/(varsigma);
z_pos = (v_pos-v_0)/(varsigma);

% in this bit is substitute the bounds in for the erf
g_neg = (a1*exp(-b1*z_neg.^2));
g_pos = ((1 - a2*exp(-b2*z_pos.^2)));

g_approx = [g_neg g_pos(2:end)];

figure
plot(v,g,v_neg,g_neg,v_pos,g_pos,v,g_approx)
