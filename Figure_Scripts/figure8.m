close all;
clear;
clc;

%{
Figure 8

The point: One can calculate the AUC ratio from fitted kinetic parameters or integrate under the time 
curves and get the same answer. This proves (a) that we did the math right in the paper and (b) that
the kinetic curves fit the data well. 
%}

sim_data = load('./Fig8Data.mat');
auc_c = sim_data.auc_c;
auc_m = sim_data.auc_m;

% Crunch the numbers for a linear fit
lin_fit = sim_data.lm; 
r2 = lin_fit.Rsquared.Ordinary;
rho = corrcoef(auc_c, auc_m);
rho = rho(1,2);
diag = 0 :0.01: max(auc_c)+.1;
fitted_line = lin_fit.Coefficients{2,1}*diag + lin_fit.Coefficients{1,1};

% Crunch the numbers for a Bland-Altman analysis
bland_x = (auc_m+auc_c) / 2;
bland_y = auc_m - auc_c;
mu = mean(bland_y, 'all');
s = std(bland_y, [], 'all');
horz_ax = [min(bland_x), max(bland_x(end))];


%% Plot
fig8 = figure('Name', 'Figure 8'); 
t = tiledlayout(1,2, 'TileSpacing','tight', 'Padding','tight');
fontsize = 15;
markersize = fontsize * 0.8;
linewidth = 1.5; 
fig8.Units = 'inches';
fig8.Position = [2.6146, 3.6354, 8.8646, 4.3750];

% Plot the linear regression
ax1 = nexttile;
axis square;
p3 = plot(auc_c, auc_m, 'kx',...
    'MarkerSize',markersize, 'LineWidth',linewidth, 'DisplayName', 'AUC Ratio');
hold on;
p2 = plot(diag, fitted_line, 'r-',...
    'LineWidth',linewidth, 'DisplayName','Linear Fit');
p1 = plot(diag, diag, 'k--',...
    'LineWidth',linewidth, 'DisplayName','Diagonal');
hold off;
xlabel('Model-Derived AUC Ratio', 'FontSize',fontsize);
ylabel('Integrated AUC Ratio', 'FontSize',fontsize);
legend([p1,p2,p3], 'FontSize',fontsize, 'Location','northwest');
ylim([0,max(diag)]);
xlim([0,max(diag)]);
% Annotations
slope_note = sprintf('y = %.2f x + %.2f', lin_fit.Coefficients{2,1}, lin_fit.Coefficients{1,1});
r2_note = sprintf('R^2 = %.2f', r2);
rho_note = sprintf('\\rho = %.2f', rho);
x_shift = 1.0;
y_shift = .15;
spacing = 0.15;
text(x_shift, 3*spacing+y_shift, slope_note, 'FontSize',fontsize, 'FontWeight','bold');
text(x_shift, 2*spacing+y_shift, r2_note, 'FontSize',fontsize, 'FontWeight','bold');
text(x_shift, spacing+y_shift,   rho_note, 'FontSize',fontsize, 'FontWeight','bold');
text(x_shift, y_shift,       'p << 0.001', 'FontSize',fontsize, 'FontWeight','bold');

% Plot the Bland-Altman style analysis
ax2 = nexttile;
axis square;
hold on;
plot(horz_ax, 1.96*s*ones(size(horz_ax)) + mu, 'r:',...
    'LineWidth',linewidth);
plot(horz_ax, -1.96*s*ones(size(horz_ax)) + mu, 'r:',...
    'LineWidth',linewidth);
plot(horz_ax, mu*ones(size(horz_ax)), 'b-',...
    'LineWidth',linewidth);
plot(bland_x, bland_y, 'kx',...
    'MarkerSize',markersize, 'LineWidth',linewidth);
hold off;
xlabel('Mean', 'FontSize',fontsize);
ylabel('Difference', 'FontSize',fontsize);
% ylim([-2.5*s+mu, 2.5*s+mu]);
box on;
grid on; 
mean_annot = sprintf('Mean: %.2f',mu);
pos_std_annot = sprintf('1.96s: %.2f',1.96*s+mu);
neg_std_annot = sprintf('-1.96s: %.2f', -1.96*s+mu);
text_offset = 0.015; 
text(max(bland_x),mu + text_offset,mean_annot,...
    'FontSize',fontsize/1.25, 'HorizontalAlignment','right');
text(max(bland_x),1.96*s+mu - text_offset,pos_std_annot,...
    'FontSize',fontsize/1.25, 'HorizontalAlignment','right');
text(max(bland_x),-1.96*s+mu + text_offset,neg_std_annot,...
    'FontSize',fontsize/1.25, 'HorizontalAlignment','right');

add_subpanel_label('A', ax1, fontsize,[-.05,-.07]);
add_subpanel_label('B', ax2, fontsize,[-.07,-.07]);


function add_subpanel_label(label, axis, fontsize, shift)
    if nargin ~= 3
        xshift = shift(1);
        yshift = shift(2);
    else
        xshift = 0;
        yshift = 0; 
    end

    pos = axis.Position;

    annotation('textbox',[pos(1)+xshift-0.01, pos(2)+pos(4)-0.03+yshift, 0.03, 0.03],...
        'String',label, 'FontSize',fontsize*2, 'FontWeight','bold',...
        'EdgeColor','none', 'VerticalAlignment','bottom');
end
