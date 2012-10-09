clc
close all
clear

FS_labels = 8;
FS_ticks = 8;

xmin = -2;
xmax = 2;
Npoints = 1000000;
x = linspace(xmin,xmax,Npoints);
% xmin:0.00001:xmax;

ymin = 0.000001;
ymax = 1-ymin;
y = linspace(ymin,ymax,Npoints);

varsigma = 1;
v_0 = 0;

% first plot the sigmoid
v = -5:0.01:5;                  % voltage range
g = 1./(1 + exp(-v));       % firing rate

figure('color','white','units','centimeters','position',[2 2 9 15],'papersize',[9 15],'filename','Distributions.pdf')

subplot(411)
plot(g,v,'k')
xlabel('Firing Rate','fontsize',FS_labels)
ylabel('Potential','fontsize',FS_labels)
box off
xlim([xmin,xmax])
set(gca,'fontsize',FS_ticks)
leg = legend('$g(v) = \left(1+e^{-v}\right)^{-1}$');
set(leg,'box','off','orientation','horizontal','interpreter','latex','position',[0.5 0.94 0 0],'fontsize',FS_labels)%'location','northoutside')

%%
% now plot the distributions for different variances
sigma_x = 0.6;

% plot the pdf of a gaussian that is the input

f_X_x = exp(-x.^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));


% find the shape of the gaussian after it has gone through the sigmoid
h_dash_y = 1./(y.*(1-y));
x = -log((1-y)./y);
f_Y_y = h_dash_y .* exp(-x.^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));

sum(f_X_x)
sum(f_Y_y)

subplot(412)
[AX,H1,H2] = plotyy(x,f_X_x,y,f_Y_y);

set(H2,'LineStyle','--')
% set(get(AX(1),'Ylabel'),'String','f(x)') 
% set(get(AX(2),'Ylabel'),'String','f(y)')
set(AX(1),'xlim',[-2 2],'ylim',[0 0.7],'xtick',[-2,0,2],'ytick',[0 0.3 0.6],'fontsize',FS_ticks)
set(AX(2),'xlim',[-2 2],'ylim',[0 3],'xtick',[-2,0,2],'ytick',[0 1 2 3],'fontsize',FS_ticks)
box off
leg = legend('$f(x) = \mathcal{N}(0,0.6^2)$','$f(y)$');
set(leg,'box','off','orientation','horizontal','interpreter','latex','position',[0.5 0.71 0 0],'fontsize',FS_labels)%'location','northoutside')

% axis tight



%%
sigma_x = 1.2;      % change the variance of the input

% plot the pdf of a gaussian that is the input

f_X_x = exp(-x.^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));    % define new input distribution

sum(f_X_x)

% find the shape of the gaussian after it has gone through the sigmoid
h_dash_y = 1./(y.*(1-y));
x = -log((1-y)./y);
f_Y_y = h_dash_y .* exp(-x.^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));

sum(f_Y_y)

subplot(413)
[AX,H1,H2] = plotyy(x,f_X_x,y,f_Y_y);

% set(get(AX(1),'Ylabel'),'String','f(x)') 
% set(get(AX(2),'Ylabel'),'String','f(y)')
set(AX(1),'xlim',[-2 2],'ylim',[0 0.7],'xtick',[-2,0,2],'ytick',[0 0.3 0.6],'fontsize',FS_ticks)
set(AX(2),'xlim',[-2 2],'ylim',[0 3],'xtick',[-2,0,2],'ytick',[0 1 2 3],'fontsize',FS_ticks)
box off
leg = legend('$f(x) = \mathcal{N}(0,1.2^2)$','$f(y)$');
set(leg,'box','off','orientation','horizontal','interpreter','latex','position',[0.5 0.48 0 0],'fontsize',FS_labels)
set(H2,'LineStyle','--')

%%
sigma_x = 1.8;      % change the variance of the input

% plot the pdf of a gaussian that is the input

f_X_x = exp(-x.^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));    % define new input distribution

sum(f_X_x)

% find the shape of the gaussian after it has gone through the sigmoid
h_dash_y = 1./(y.*(1-y));
x = -log((1-y)./y);
f_Y_y = h_dash_y .* exp(-x.^2 / (2*sigma_x^2)) / (sigma_x * sqrt(2*pi));

sum(f_Y_y)

subplot(414)
[AX,H1,H2] = plotyy(x,f_X_x,y,f_Y_y);

% set(get(AX(1),'Ylabel'),'String','f(x)') 
% set(get(AX(2),'Ylabel'),'String','f(y)')
set(AX(1),'xlim',[-2 2],'ylim',[0 0.7],'xtick',[-2,0,2],'ytick',[0 0.3 0.6],'fontsize',FS_ticks)
set(AX(2),'xlim',[-2 2],'ylim',[0 3],'xtick',[-2,0,2],'ytick',[0 1 2 3],'fontsize',FS_ticks)
box off
set(H2,'LineStyle','--')
leg = legend('$f(x) = \mathcal{N}(0,1.8^2)$','$f(y)$');
set(leg,'box','off','orientation','horizontal','interpreter','latex','position',[0.5 0.26 0 0],'fontsize',FS_labels)