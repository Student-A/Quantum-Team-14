function plots()

fs = 20;
fs2 = 20;
data = load('E_f');
k = size(data);
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
title('Upper Bound for E_f','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'E_f.eps','psc2');

 s3 = ['Probability_Distribution'];
 data = load(s3);
 
 [X Y] = meshgrid([floor(min(data(:,1))):1: ceil(max(data(:,1)))],[floor(min(data(:,2))):1: ceil(max(data(:,2)))]);
 
  Z = griddata(data(:,1), data(:,2), data(:,3),X,Y);
    h =  figure
    surfc(X,Y,Z,'FaceColor', 'interp');
    xlabel('x','fontsize',fs);
    ylabel('y','fontsize',fs);
    %zlabel(['p_{12}(i,j;' num2str(k) ')'],'fontsize',fs);
    title('p_{12}(i,j;n)','fontsize',fs);
    colormap(hsv);
    axis([min(data(:,1)) max(data(:,1)) min(data(:,2)) max(data(:,2)) 0 max(data(:,3))])
    set(gca,'fontsize',fs2)
    saveas(h,'probability_plane.eps','psc2');
    
    
s3 = ['Average_Probability_Distribution'];
data = load(s3);
[X Y] = meshgrid([floor(min(data(:,1))):1: ceil(max(data(:,1)))],[floor(min(data(:,2))):1: ceil(max(data(:,2)))]);
 Z = griddata(data(:,1), data(:,2), data(:,3),X,Y);
    h = figure
    surfc(X,Y,Z,'FaceColor', 'interp');
    xlabel('x','fontsize',fs);
    ylabel('y','fontsize',fs);
    %zlabel('\Pi_{12}(i,j)','fontsize',fs);
    title('\pi_{12}(i,j)','fontsize',fs);
    colormap(hsv);
      axis([min(data(:,1)) max(data(:,1)) min(data(:,2)) max(data(:,2)) 0 max(data(:,3))])
    set(gca,'fontsize',fs2)
    saveas(h,[s3 '.eps'],'psc2');



data = load('H_x');
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('H(x_1)','fontsize',fs);
title('H(x_1)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'H_x.eps','psc2');


data = load('H_y');
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('H(x_2)','fontsize',fs);
title('H(x_2)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'H_y.eps','psc2');


data = load('H_xy')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('H(x_1,x_2)','fontsize',fs);
title('H(x_1,x_2)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'H_xy.eps','psc2');


data = load('I_xy');
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('I(x_1;x_2)','fontsize',fs);
title('I(x_1;x_2)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'I_xy.eps','psc2');


data = load('S_x');
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('S(\rho_{P,1})','fontsize',fs);
title('S(\rho_{P,1})','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'S_x.eps','psc2');



data = load('S_y')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('S(\rho_{P,2})','fontsize',fs);
title('S(\rho_{P,2})','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'S_y.eps','psc2');


data = load('Von_Newman_entropy of coin state')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('E_{C,P}','fontsize',fs);
title('E_{C,P}','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'S_C.eps','psc2');


data = load('Iv_xy')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('I(\rho_{P,12})','fontsize',fs);
title('I(\rho_{P,12})','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'Iv_xy.eps','psc2');




data = load('Quantum Discord')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('\delta(x_2:x_1)','fontsize',fs);
title('\delta(x_2:x_1)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'QD.eps','psc2');



data = load('xy_cov')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('cov_{1,2}(x_1,x_2)','fontsize',fs);
title('cov(x_1,x_2)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'xy_cov.eps','psc2');



data = load('Average hitting time')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('N_a^{(1)}(i_0)','fontsize',fs);
title('N_a^{(1)}(i_0)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'av_hit.eps','psc2');

data = load('one Shot probability to hit')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('P_o^1(i_0;n)','fontsize',fs);
title('P_o^1(i_0;n)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'one_shot_hit.eps','psc2');

data = load('Concurrence hitting time')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('P_c^{(1)}(i_0;n)','fontsize',fs);
title('P_c^{(1)}(i_0;n)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'concurrent.eps','psc2');




data = load('First time to hit')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('P_f^{(1)}(i_0,n)','fontsize',fs);
title('P_f^{(1)}(i_0,n)','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'first_time_hit.eps','psc2');





data = load('mean_dist')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('<x-y>','fontsize',fs);
title('<|X_1-X_2|>','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'mean_dist.eps','psc2');


data = load('mean_x')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('<x>','fontsize',fs);
title('<x_1>','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'mean_x.eps','psc2');

data = load('mean_y')
h = figure
plot(data(:,1), data(:,2));
xlabel('time','fontsize',fs);
%ylabel('<y>','fontsize',fs);
title('<x_2>','fontsize',fs);
set(gca,'fontsize',fs2)
saveas(h,'mean_y.eps','psc2');
