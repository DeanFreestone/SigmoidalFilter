% test the square of two error function

clc
clear
close all

min_w = 0;
dw = 0.001;
max_w = 20;

w = min_w:dw:max_w;


dx = 0.001;
x = -10:dx:10;
errorfunction = zeros(1,length(x));

for n=1:length(x)    
    errorfunction(n) = sum(exp(-(w-x(n)).^2/2)/sqrt(2*pi)) * dw;
end

figure
plot(x,errorfunction,x,0.5*(erf(x/sqrt(2))+1))
legend('numerically int erf','matlab erf')

figure
plot(x,errorfunction.^2,x,(0.5*(erf(x/sqrt(2))+1)).^2)
legend('squared numerically int erf','squared matlab erf')