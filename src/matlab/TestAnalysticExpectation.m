% clear
% clc
% close all
% 
% mu = 0;
% sigma = .5;
% x = -10:0.001:10;
% cummdensfunc = (erf((1/sigma)*(x-mu) / sqrt(2)) + 1) / 2;
% 
% % g_x2 = 1./(1+exp(sigma*(mu-x)));
% % 
% % figure,plot(x,cummdensfunc,x,g_x2)
% 
% % define the gaussian and the realisations of the RV
% N_realisations = 1e7;
% mu_x = 0;
% sigma_x = 2;
% x = mu_x + sigma_x*randn(1,N_realisations);
% 
% % map through the sigmoid
% g_x = (erf((1/sigma)*(x-mu) / sqrt(2)) + 1) / 2;
% % g_x2 = 1./(1+exp(varsigma*(mu-x)));
%     
% % figure
% % plot(g_x,'.k')
% mean(g_x)
% 
% y = mu_x/sqrt(1+sigma_x^2); % use the stats o f the data
% mu_analytic = (erf((1/sigma)*(y-mu) / sqrt(2)) + 1) / 2 % pass this through sigmoid

% this is from the paper, demonstrating that it works for a simple sigmoid
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear
clc
close all

FS_Label = 8;
FS_ticks = 8;

mu = 1;
s = 1;

N_realisations = 1e8;
x = mu + s*randn(1,N_realisations);
y = 0.5*(erf(x/sqrt(2)) + 1);           % transform RV through error function
numerical_M = mean(y);
numerical_V = var(y);
disp(['Numerical Mean = ' num2str(numerical_M)])
disp(['Numerical Variance = ' num2str(numerical_V)])


z = mu/sqrt(1+s^2);
M = 0.5*(erf(z/sqrt(2)) + 1);
disp(['Analytic Mean = ' num2str(M)])

V1 = mean(y.^2) - M^2;
disp(['E(y^2) - M^2, Variance = ' num2str(V1)])

V_temp = mu/sqrt(1+s^2);
Phi_V_temp = 0.5*(erf(V_temp/sqrt(2)) + 1);
V = Phi_V_temp*(1-Phi_V_temp)*s^2/(s^2 + pi/2);
disp(['Analytic Variance = ' num2str(V)])

E_y2 = erfc(mu/sqrt(4*s^2+2))^2*(2*s^2 + 1)/(4*sqrt(2*s^2+1));
% 
V2 = 1-2*M + E_y2 - M^2;
disp(['Analytic Variance = ' num2str(V2)])

%%
% NBINS = 50;
% 
% figure('color','white','units','centimeters','position',[2 2 9 4],'papersize',[9 4],'filename','Distributions.pdf')
% 
% AX(1) = subplot(121);
% [Nx Xx] = hist(x,NBINS);
% bar(Xx,Nx/max(Nx))
% xlabel('Bin Center','fontsize',FS_labels)
% ylabel('Norm Bin Count','fontsize',FS_labels)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% set(gca,'fontsize',FS_ticks)
% 
% xlim([-2 4])
% box off
% %%
% AX(2) = subplot(122);
% [Ny Xy] = hist(y,NBINS);
% bar(Xy,Ny/max(Ny))
% set(gca,'fontsize',FS_ticks)
% 
% xlabel('Bin Center','fontsize',FS_labels)
% ylabel('Norm Bin Count','fontsize',FS_labels)
% xlim([-0 1])
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','w','EdgeColor','k')
% box off
% hold
% plot([M M],[0 1],'r','linewidth',1.5)
% plot([numerical_M numerical_M],[0 1],'b','linewidth',0.5)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%
% now if we make the sigmoid like a NMM see what happens
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mu_s = 6;
s_s = 3.033898556455646;

% pass the realisations of the RV through sigmoid
y = 0.5*(erf((x-mu_s)/(s_s*sqrt(2))) +1);
mean_y = mean(y);
disp(['The Numerical Mean of the NMM Like Error Sigmoid is: ' num2str(mean_y)])
z = (mu-mu_s)/sqrt(s_s^2+s^2);
M = 0.5*(erf(z/sqrt(2)) + 1);
disp(['The Analytic Mean of the NMM Like Error Sigmoid is: ' num2str(M)])

%%
%plot the sigmoidal function
v_0 = 6;            % threshold
vmin = v_0-15;
vmax = v_0+15;
v = linspace(vmin,vmax,1000);
varsigma_l = 0.56;       % slope
NIterations = 10000;
RMSE = zeros(1,NIterations);
varsigma_e = linspace(3.033898,3.033899,NIterations);
for n=1:NIterations    
    g_e = 0.5*(erf((v-v_0)/(varsigma_e(n)*sqrt(2))) +1);
    g_l = 1 ./ (1 + exp(varsigma_l*(v_0 - v)));
    error = g_l-g_e;
    RMSE(n) = sqrt(mean(error.^2));
end
figure,plot(varsigma_e,RMSE)
axis tight
[y I] = min(RMSE);
varsigma_e = varsigma_e(I);
g_e = 0.5*(erf((v-v_0)/(varsigma_e*sqrt(2))) +1);
g_l = 1 ./ (1 + exp(varsigma_l*(v_0 - v)));

figure('color','white','units','centimeters','position',[2 2 9 7],'papersize',[7 9],'filename','Mapping.pdf')
subplot(311)
plot(v,g_l,'k')
xlim([vmin,vmax])
ylim([-0.01 1.01])
box off
set(gca,'fontsize',FS_ticks)
ylabel('$g_l(v)$','fontsize',FS_Label,'interpreter','latex')

subplot(312)
plot(v,g_e,'k')
xlim([vmin,vmax])
ylim([-0.01 1.01])
box off
set(gca,'fontsize',FS_ticks)
ylabel('$g_e(v)$','fontsize',FS_Label,'interpreter','latex')

subplot(313)
error = g_l-g_e;
plot(v,error,'k')
xlim([vmin,vmax])
ylim([-0.01 0.01])
box off
set(gca,'fontsize',FS_ticks)
xlabel('Mean Membrane Potential (mV)','fontsize',FS_Label)
ylabel('$g_l(v) - g_e(v)$','fontsize',FS_Label,'interpreter','latex')