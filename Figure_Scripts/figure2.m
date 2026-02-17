close all;
clear; 
clc; 

addpath('../Simulation_Functions');
%{
Figure 2

The point: The new model (3PCs or P2L3s) can fit time curves from the full
model (3PC or P2L3) with very little error. That means 3PCs is a good
approximation of the full system.
%}


%% Generate noiseless time curves from three physiological regimes
% Acquisition parameters
fdv.ntp = 64;                          % Number of time points
fdv.NSeg = 1;                          % Segments per time point
fdv.NFlips = fdv.ntp * fdv.NSeg;       % Total number of excitations
fdv.TE = 0; %0.0219;                   % bbEPI; set TE=0 to ignore T2* differences
fdv.TR = 3 * ones(1,fdv.NFlips);
fdv.taxis = cumsum(fdv.TR) - fdv.TR(1); 
fdv.FlipAngle = repmat( [20;30;20;30], [1,fdv.NFlips] ); 
fdv.verbose = 0; 
fdv.SPBins = 1;                        % Description of slice profiles
fdv.SPWeights = 1;                     % Can be turned to 1 for block pulse

% Physics parameters
base_knowns    = {'T1Pyr', 'T1Lac', 'R1Lx', 'T2spv', 'T2spe', 'T2spc', 'T2slc', 'Pe0', 'Pi0', 'Le0', 'Li0', 'FAp', 'FAl'};
base_knownvals = [  43,       33,    0.0,     0.1,    0.055,   0.055,   0.033,    0,     0,     0,     0,     20,   30];
% VIF Parameters
base_knowns    = {base_knowns{:}, 'Gam1', 'Gam2', 'tdel'};
base_knownvals = [base_knownvals,  2.8,     4.5,    10];
fdv.VIFP = gampdf(fdv.taxis-6, 2.8, 4.5); 
fdv.VIFL = zeros( 1,length(fdv.VIFP) ); 
fdv.UseVIF = 1;
fprintf('\nPyruvate locales for Full 3PC model...\n');

% Physiological parameters for dominant Pyr IV
fdv.knowns    = {base_knowns{:}, 'vb', 'vef', 'klp', 'kecl'};
fdv.knownvals = [base_knownvals, 0.15,  0.5,  0.02,   0.03];
fdv.fitvars   = {'kpl', 'VIFScale', 'kve', 'kecp', 'R1Lx'};
iv_parms      = [ 0.5,     1000,    0.011,  0.03,     0 ]; 
% Lowering extravasation and increasing vascular volume means more Pyr IV
% Generate synthetic data
[EV, IV, vols] = P2L3(iv_parms, fdv);
iv_data = [
    IV(1,:)*vols(1) + EV(1,:)*vols(2) + EV(3,:)*vols(3);  % Pyruvate signal
    EV(2,:)*vols(2) + EV(4,:)*vols(3) ];                  % Lactate signal
% Confirm that we get a dominant signal from EE
pvfsv = P2L3_comp_contr(iv_parms, fdv);
fprintf('IV Dominant: Pyr IV = %4.1f%% | EE = %4.1f%% |  C = %4.1f%% \n', pvfsv*100);

% Physiological parameters for dominant Pyr EE
fdv.knowns    = {base_knowns{:}, 'vb', 'vef', 'klp', 'kecl'};
fdv.knownvals = [base_knownvals,  0.1,  0.5,  0.02,   0.03 ];
%fdv.fitvars= {'kpl', 'VIFScale', 'kve', 'kecp', 'R1Lx'};
ee_parms    = [ 0.5,     1000,    0.066,  0.03,   0.0];
% Generate synthetic data
[EV, IV, vols] = P2L3(ee_parms, fdv);
ee_data = [
    IV(1,:)*vols(1) + EV(1,:)*vols(2) + EV(3,:)*vols(3);  % Pyruvate signal
    EV(2,:)*vols(2) + EV(4,:)*vols(3) ];                  % Lactate signal
% Confirm that we get a dominant signal from EE
pvfse = P2L3_comp_contr(ee_parms, fdv);
fprintf('EE Dominant: Pyr IV = %4.1f%% | EE = %4.1f%% |  C = %4.1f%% \n', pvfse*100);

% Physiological parameters for dominant Pyr C
fdv.knowns    = {base_knowns{:}, 'vb', 'vef', 'klp', 'kecl'};
fdv.knownvals = [base_knownvals, 0.02,  0.2, 0.0008,   0.2 ];
% fdv.fitvars = {'kpl', 'VIFScale', 'kve', 'kecp', 'R1Lx'};
c_parms      = [  0.2,    1000,      0.08,  0.2,     0.0 ]; 
% Generate synthetic data
[EV, IV, vols] = P2L3(c_parms, fdv);
c_data = [
    IV(1,:)*vols(1) + EV(1,:)*vols(2) + EV(3,:)*vols(3);  % Pyruvate signal
    EV(2,:)*vols(2) + EV(4,:)*vols(3) ];                  % Lactate signal
% Confirm that we get a dominant signal from EE
pvfsc = P2L3_comp_contr(c_parms, fdv);
fprintf(' C Dominant: Pyr IV = %4.1f%% | EE = %4.1f%% |  C = %4.1f%% \n\n', pvfsc*100);


%% Fit the 3PCs model to all three regimes
n_fits = 500;
verbose_fitting = 0;
base_knowns    = {'T1Pyr', 'T1Lac', 'klp', 'kecl', 'T2spv', 'T2spe', 'T2spc', 'T2slc', 'Pe0', 'Pi0', 'Le0', 'Li0', 'FAp', 'FAl'};
base_knownvals = [  43,       33,    0.0,   0.0,     0.1,    0.055,   0.055,   0.033,    0,     0,     0,     0,     20,   30];
fdv.FlipAngle = [20,0; 0,30] * ones(2,fdv.NFlips);
fprintf('\nPyruvate locales for fitted 3PCs curves...\n');

% Fit Pyr IV dominant
fdv.knowns    = {base_knowns{:}, 'vb', 'vef'}; 
fdv.knownvals = [base_knownvals, 0.15,  0.5];
fdv.data = iv_data; 
[iv_popt, iv_fitted, fdv_modded, iv_nmse, iv_hits] = fit_3PCs(fdv, n_fits, verbose_fitting);
[~,iv_pyr_fracs] = P2L3sAUC(iv_popt, fdv_modded,fdv_modded);
fprintf('IV Dominant: Pyr IV = %4.1f%% | EE = %4.1f%% |  C = %4.1f%% \n', iv_pyr_fracs*100);
% Fit Pyr EE dominant
fdv.knowns    = {base_knowns{:}, 'vb', 'vef'}; 
fdv.knownvals = [base_knownvals,  0.1,  0.5];
fdv.data = ee_data;
[ee_popt, ee_fitted, fdv_modded, ee_nmse, ee_hits] = fit_3PCs(fdv, n_fits, verbose_fitting);
[~,ee_pyr_fracs] = P2L3sAUC(ee_popt, fdv_modded,fdv_modded);
fprintf('EE Dominant: Pyr IV = %4.1f%% | EE = %4.1f%% |  C = %4.1f%% \n', ee_pyr_fracs*100);
% Fit Pyr C dominant
fdv.knowns    = {base_knowns{:}, 'vb', 'vef'}; 
fdv.knownvals = [base_knownvals, 0.02,  0.2];
fdv.data = c_data;
[c_popt, c_fitted, fdv_modded, c_nmse, c_hits] = fit_3PCs(fdv, n_fits, verbose_fitting); 
[~,iv_pyr_fracs] = P2L3sAUC( c_popt, fdv_modded,fdv_modded);
fprintf(' C Dominant: Pyr IV = %4.1f%% | EE = %4.1f%% |  C = %4.1f%% \n\n', iv_pyr_fracs*100);


%% Format a pretty picture

% Normalize to tallest synthetic pyruvate, just for display
pyr_peaks = [
    max(iv_data, [], 'all');
    max(ee_data, [], 'all'); 
    max( c_data, [], 'all')];
norm_factor = max(pyr_peaks); 

iv_data = iv_data / norm_factor;
ee_data = ee_data / norm_factor;
 c_data =  c_data / norm_factor; 

iv_fitted = iv_fitted / norm_factor;
ee_fitted = ee_fitted / norm_factor;
 c_fitted =  c_fitted / norm_factor; 

% Plotting
fig2 = figure('Name', 'Figure2');
tiledlayout(1,3, 'TileSpacing', 'compact');
fontsize = 15;
linewidth = 1.25;
markersize = 12;

fig2.Units = 'inches';
fig2.Position = [0.5, 5.5, 6.9 *2, 4.375];

% Pyr IV dominates
plot_y_max = max(iv_data, [], 'all');
ax1 = nexttile;
plot(fdv.taxis, iv_data(1,:), 'gx', 'DisplayName', '3PC Simulated Pyruvate',...
    'MarkerSize', markersize);
hold on;
plot(fdv.taxis, iv_data(2,:), 'b.', 'DisplayName', '3PC Simulated Lactate',...
    'MarkerSize', markersize*1.15);
plot(fdv.taxis, iv_fitted(1,:), 'g-', 'DisplayName', '3PCs Fitted Pyruvate',...
    'MarkerSize', markersize, 'LineWidth', linewidth);
plot(fdv.taxis, iv_fitted(2,:), 'b--', 'DisplayName', '3PCs Fitted Lactate', ...
    'MarkerSize', markersize, 'LineWidth', linewidth);
hold off; 
grid on;
add_subpanel_label('A', ax1, fontsize, [0,-0.05]);
xlabel('Time (s)', 'FontSize', fontsize); 
ylabel('Normalized Synthetic Signal', 'FontSize', fontsize);
title('Pyr_{iv} Dominant', 'FontSize', fontsize);
legend('FontSize', fontsize*0.8);
ylim([0, plot_y_max]);
xlim([-1, 150]);
nmse_label = sprintf('NMSE: %.2e', iv_nmse);
text(60, plot_y_max*0.5, nmse_label, 'FontSize', fontsize*1.1)

% Pyr EE dominates
plot_y_max = max(ee_data, [], 'all');
ax2 = nexttile;
plot(fdv.taxis, ee_data(1,:), 'gx',...
    'MarkerSize', markersize);
hold on;
plot(fdv.taxis, ee_data(2,:), 'b.', 'DisplayName', '3PC Simulated Lactate',...
    'MarkerSize', markersize*1.15);
plot(fdv.taxis, ee_fitted(1,:), 'g-', 'DisplayName', '3PCs Fitted Pyruvate',...
    'MarkerSize', markersize, 'LineWidth', linewidth);
plot(fdv.taxis, ee_fitted(2,:), 'b--', 'DisplayName', '3PCs Fitted Lactate', ...
    'MarkerSize', markersize, 'LineWidth', linewidth);
hold off; 
grid on;
add_subpanel_label('B', ax2, fontsize);
xlabel('Time (s)', 'FontSize', fontsize); 
title('Pyr_{ee} Dominant', 'FontSize', fontsize);
ylim([0, plot_y_max]);
xlim([-1, 150]);
nmse_label = sprintf('NMSE: %.2e', ee_nmse);
text(60, plot_y_max*0.5, nmse_label, 'FontSize', fontsize*1.1)

% Pyr C dominates
plot_y_max = max(c_data, [], 'all'); 
ax3 = nexttile;
plot(fdv.taxis, c_data(1,:), 'gx',...
    'MarkerSize', markersize);
hold on;
plot(fdv.taxis, c_data(2,:), 'b.',...
    'MarkerSize', markersize*1.15);
plot(fdv.taxis, c_fitted(1,:), 'g-',...
    'MarkerSize', markersize, 'LineWidth', linewidth);
plot(fdv.taxis, c_fitted(2,:), 'b--',...
    'MarkerSize', markersize, 'LineWidth', linewidth);
hold off; 
grid on;
add_subpanel_label('C', ax3, fontsize);
xlabel('Time (s)', 'FontSize', fontsize); 
title('Pyr_c Dominant', 'FontSize', fontsize);
ylim([0, plot_y_max]);
xlim([-1, 150]);
nmse_label = sprintf('NMSE: %.2e', c_nmse);
text(60, plot_y_max*0.5, nmse_label, 'FontSize', fontsize*1.1)

axes = [ax1, ax2, ax3];
for i = 1:length(axes)
    axes(i).GridColor = 'black';
    axes(i).FontSize = fontsize / 1.5;
    axes(i).XLabel.FontSize = fontsize;
    axes(i).YLabel.FontSize = fontsize;
    axes(i).Title.FontSize = fontsize * 1.25;
end

% Export as MATLAB fig file for later editing
savefig(fig2, '../Media/figure2.fig');

% Export as vector image for publication
exportgraphics(fig2, '../Media/figure2.eps',...
    'Resolution', 300,...
    'BackgroundColor', 'white',...
    'ContentType', 'vector');

% Export as pixel image for preview
exportgraphics(fig2, '../Media/figure2.png',...
    'Resolution', 300,...
    'BackgroundColor', 'white');


%% Display fit convergence to make sure that we're using enough guesses
% Plot the number of new hits over the course of fitting
figure('Name','Convergence Study pt1', 'Position',[1377, 440, 509, 379]);
tiledlayout(2,1, 'TileSpacing','tight');

iv_conv = [0, diff( iv_hits(1,:) )];
ee_conv = [0, diff( ee_hits(1,:) )];
 c_conv = [0, diff(  c_hits(1,:) )];

iv_resid = iv_hits(2,:);
ee_resid = ee_hits(2,:);
 c_resid =  c_hits(2,:);
 
loglog(iv_hits(1,:), iv_conv, 'x-',...
    'color',[0, 0.4470, 0.7410], 'DisplayName','Pyr_{iv} Dominant');
hold on;
loglog(ee_hits(1,:), ee_conv, 'o-',...
    'color',[0.6350 0.0780 0.1840], 'DisplayName','Pyr_{ee} Dominant');
loglog( c_hits(1,:),  c_conv, 'o-', 'MarkerFaceColor',[0.9290 0.6940 0.1250],...
    'color',[0.9290 0.6940 0.1250], 'DisplayName','Pyr_c Dominant');
hold off; 
legend('FontSize',fontsize*0.8, 'Location','southeast');
grid on;
ylabel('Guesses since last hit', 'FontSize',fontsize);
xlim([0,n_fits]);
xlabel('Nth Guess', 'FontSize',fontsize);

% Plot the residual over the course of fitting
figure('Name','Convergence Study pt2', 'Position',[48, 66, 1324, 360]);
tiledlayout(1,3, 'TileSpacing','tight');

nexttile;
plot(iv_hits(1,:), iv_resid, 'x-',...
    'color',[0, 0.4470, 0.7410], 'DisplayName','Pyr_{iv} Dominant');
grid on;
xlabel('Nth Guess', 'FontSize',fontsize);
ylabel('Residual Error', 'FontSize',fontsize);

nexttile;
plot(ee_hits(1,:), ee_resid, 'o-',...
    'color',[0.6350 0.0780 0.1840], 'DisplayName','Pyr_{ee} Dominant');
xlabel('Nth Guess', 'FontSize',fontsize);
grid on; 

nexttile;
plot( c_hits(1,:),  c_resid, 'o-', 'MarkerFaceColor',[0.9290 0.6940 0.1250],...
    'color',[0.9290 0.6940 0.1250], 'DisplayName','Pyr_c Dominant');
grid on;
xlabel('Nth Guess', 'FontSize',fontsize);


%% Local function for fitting
function [popt, fit_curves, fdv_modded, nmse, convg_study] = fit_3PCs(fdv, n_fits, verbose)
    jbopts = optimset('display', 'off'); 
    best_resid = Inf;
    best_fits = 0;
    
    if verbose
        fprintf('Starting fitting with %d repetitions...', n_fits);
        tic;
    end

    % % For fitting T1s without extra loss factor
    % fdv.fitvars   = {'kpl', 'VIFScale', 'kve', 'kecp', 'vb', 'vef', 'T1Pyr', 'T1Lac'};
    % fdv.knowns    = {'klp', 'kecl', 'R1Lx', 'T2spv', 'T2spe', 'T2spc', 'T2slc', 'Pe0', 'Pi0', 'Le0', 'Li0', 'FAp', 'FAl'};
    % fdv.knownvals = [ 0.0,   0.0,    0.0,     0.1,    0.055,   0.055,   0.033,    0,     0,     0,     0,     20,   30];
    % UB            = [ 10,      Inf,       10,    10,    1,     1,     Inf,      Inf ];
    % LB            = [  0,       0,         0,     0,    0,     0,      0,        0 ];
    % p0            = [ 0.5,     1000,     0.07,  0.02,   0.1,  0.5,    43,       33 ];
    
    % % For fitting extra loss factor without T1s
    fdv.fitvars   = {'kpl', 'VIFScale', 'kve', 'kecp', 'R1Lx'};
    UB            = [  1,      Inf,       1,     1,       1 ];
    LB            = [  0,       0,         0,    0,       0 ];
    p0            = [ 0.5,     1000,     0.07,  0.02,     0 ];

    convergence = 0;
    resids = Inf; 
    for ii=1:n_fits
        guess = (p0*9).*rand(1,length(p0));
        [fits, resid] = lsqnonlin(@(x) P2L3sErr(x, fdv), guess, LB, UB, jbopts);

        if resid < best_resid
            best_resid = resid;
            best_fits = fits;
            convergence = [convergence, ii];
            resids = [resids, best_resid];
        end
    end
    convg_study = [convergence; resids];

    popt = best_fits;
    [EV, IV, vols] = P2L3s(popt, fdv);
    fit_curves = [
        IV*vols(1) + EV(1,:)*vols(2) + EV(2,:)*vols(3);  % Pyruvate signal
        EV(3,:)*vols(3)                                  % Lactate signal
    ]; 
    fdv_modded = fdv;
    nmse = calc_nmse(fdv.data, fit_curves);

    if verbose
        fprintf('finished. Elapsed time: %.2f min\n', toc()/60);
    end
end

% Local function for adding subpanel labels
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

