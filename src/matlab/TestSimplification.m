clear
clc
close all

mu = 1;
s = 2;

dw = 0.001;
dz = dw;

w = 0:0.001:30;
z = w;

temp = zeros(1,length(z));

coeff = (s^2+1)/(4*s^2+2);

for n=1:length(z)
    temp(n) = sum(exp(-coeff*(w - (mu+s^2*z(n))/(s^2+1)).^2))*dw;
end

temp2 = sqrt(pi*(1-1/(2*s^2+2))) - sqrt(pi/2)*sqrt(2*s^4+3*s^2+1)/(s^2+1)*erf((-z*s^2-mu)/sqrt(4*s^4+6*s^2+2));

figure
plot(z,temp,z,temp2)
legend('num int','analytic int')

%%
temp = zeros(1,length(z));
temp2 = zeros(1,length(z));
for n=1:length(z)
    temp(n) = sum(exp((2*mu*w - (s^2+1)*w.^2 + 2*s^2*z(n)*w)/(4*s^2+2)))*dw;
    
    temp2(n) = sum(exp(-(s^2+1)/(4*s^2+2) * (w - (mu + s^2*z(n))/(s^2+1) ).^2 ))*dw;
end

coeff = exp((mu+s^2*z).^2/(2*(s^2+1)*(2*s^2+1)));

figure
plot(z,temp,z,coeff.*temp2)
legend('before completing square','after completing square')