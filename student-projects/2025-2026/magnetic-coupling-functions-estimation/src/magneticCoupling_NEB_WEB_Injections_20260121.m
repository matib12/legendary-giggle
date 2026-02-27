run  /users/swinkels/deploy/MatlabVirgoTools/trunk/startup.m
format long g

ENV_CF_repo_path='/users/lunghini/ENV_repo/couplingfunctions/O4/Magnetic/CF_files';

FFTlength   = 300; 
overlapPerc = 50;
percToPlot  = [10,50,90];
colorForPercentiles=["red";"green";"blue"];

saveplots  = false;
resolution = 600;

%% Load other CFs
load("CF_lines_lowfreq_NEB-WEB-CEB_old.mat")
load("WEB_CF_Lines_20251211.mat")
NEB_percentiles=load_CF_percentiles_from_ENV_git(ENV_CF_repo_path,'magnetic','NEB',{'20250430','20250515'},percToPlot);
CEB_percentiles=load_CF_percentiles_from_ENV_git(ENV_CF_repo_path,'magnetic','CEB',{'20250430','20250515'},percToPlot);
WEB_percentiles=load_CF_percentiles_from_ENV_git(ENV_CF_repo_path,'magnetic','WEB',{'20250430','20250515'},percToPlot);

%%
% NEB_Metatron_outfile_old = "/virgoData/NoiseInjections/MagneticInjectionsO4/output/MagneticLine_NOISE_MAG_NEB-1453073344.txt";
NEB_Metatron_outfile = "/virgoData/NoiseInjections/MagneticInjectionsO4/output/MagneticLine_NOISE_MAG_NEB-1453290488.txt";
% CEB_Metatron_outfile = "/virgoData/NoiseInjections/MagneticInjectionsO4/output/MagneticLine_NOISE_MAG_CEB-1453302290.txt";
CEB_Metatron_outfile = "/virgoData/NoiseInjections/MagneticInjectionsO4/output/MagneticLine_NOISE_MAG_CEB-1455655608.txt";


if ~exist('NEB','var');NEB = struct();end
if ~exist('CEB','var');CEB = struct();end


% NEB=compute_CForUL_from_Metatron_outfile(NEB,NEB_Metatron_outfile,FFTlength,overlapPerc,1.8,0);
CEB=compute_CForUL_from_Metatron_outfile(CEB,CEB_Metatron_outfile,FFTlength,overlapPerc,1.8,0,CEB_percentiles);
return

%% Plot NEB Lines new + O3
set_plot_default_properties();
% close all
figure()
hold on 
plot(CF_lines_lowfreq_Old.NEB.FREQHz,CF_lines_lowfreq_Old.NEB.CFmT,'or','MarkerFaceColor','r','DisplayName','O3 Lines Injection')
plot(NEB.DATA.RDS.CForUL.freqCF,NEB.DATA.RDS.CForUL.CF,'sk','MarkerFaceColor','k','DisplayName','Lines 20260124')
plot(NEB.DATA.RDS.CForUL.freqUL,NEB.DATA.RDS.CForUL.UL, 'dg', 'MarkerFaceColor','g','DisplayName','UL Lines 20260124')
xlabel("Frequency [Hz]")
ylabel("CF [m/T]")
title("NEB Magnetic Coupling")
adjustAxisLimits(0.05,'X')
adjustAxisLimits(0.05,'Y')
ylim([1e-10 5e-4])
legend()

if saveplots
    exportgraphics(gcf, sprintf("MAG_CF_plots_Metatron/NEB_CF_Magnetic_Lines_20260124.png"), 'Resolution', resolution);
end

%% Plot CEB Lines new + O3
set_plot_default_properties();
% close all
figure()
hold on 
plot(CF_lines_lowfreq_Old.CEB.FREQHz,CF_lines_lowfreq_Old.CEB.CFmT,'or','MarkerFaceColor','r','DisplayName','O3 Lines Injection')
plot(CEB.DATA.RDS.CForUL.freqCF,CEB.DATA.RDS.CForUL.CF,'sk','MarkerFaceColor','k','DisplayName','CF Lines 20260220')
plot(CEB.DATA.RDS.CForUL.freqUL,CEB.DATA.RDS.CForUL.UL, 'dg', 'MarkerFaceColor','g','DisplayName','UL Lines 20260220')
xlabel("Frequency [Hz]")
ylabel("CF [m/T]")
title("CEB Magnetic Coupling")
adjustAxisLimits(0.05,'X')
adjustAxisLimits(0.05,'Y')
ylim([1e-10 5e-4])
legend()
if saveplots
    exportgraphics(gcf, sprintf("MAG_CF_plots_Metatron/CEB_CF_Magnetic_Lines_20260124.png"), 'Resolution', resolution);
end

%% Plot WEB Lines new + O3
set_plot_default_properties();
% close all
figure()
hold on 
plot(CF_lines_lowfreq_Old.WEB.FREQHz,CF_lines_lowfreq_Old.WEB.CFmT,'or','MarkerFaceColor','r','DisplayName','O3 Lines Injection')
plot(WEB_CForUL_20251211.freq,WEB_CForUL_20251211.CF,'sk','MarkerFaceColor','k','DisplayName','Lines 20251211')
xlabel("Frequency [Hz]")
ylabel("CF [m/T]")
title("WEB Magnetic Coupling")
adjustAxisLimits(0.05,'X')
adjustAxisLimits(0.05,'Y')
ylim([1e-10 5e-4])
legend()
if saveplots
    exportgraphics(gcf, sprintf("MAG_CF_plots_Metatron/WEB_CF_Magnetic_Lines_20251211.png"), 'Resolution', resolution);
end
% return
%% Plot NEB Lines new + all Sweeps
set_plot_default_properties();
% close all
figure()
hold on
for i=1:numel(colorForPercentiles)
    plot(NEB_percentiles.freq, NEB_percentiles.percCF(i,:),'.','Color',colorForPercentiles(i),'DisplayName',sprintf("%d %% Sweep percentile",percToPlot(i)))
end
plot(NEB.DATA.RDS.CForUL.freqCF,NEB.DATA.RDS.CForUL.CF,'sk','MarkerFaceColor','k','DisplayName','Lines 20260124')
xlabel("Frequency [Hz]")
ylabel("CF [m/T]")
title("NEB Magnetic Coupling")
adjustAxisLimits(0.05,'X')
adjustAxisLimits(0.05,'Y')
legend()
if saveplots
    exportgraphics(gcf, sprintf("MAG_CF_plots_Metatron/NEB_CF_Magnetic_Lines_20260124_plus_Sweeps.png"), 'Resolution', resolution);
end

%% Plot CEB Lines new + all Sweeps
set_plot_default_properties();
% close all
figure()
hold on
for i=1:numel(colorForPercentiles)
    plot(CEB_percentiles.freq, CEB_percentiles.percCF(i,:),'.','Color',colorForPercentiles(i),'DisplayName',sprintf("%d %% Sweep percentile",percToPlot(i)))
end
plot(CEB.DATA.RDS.CForUL.freqCF,CEB.DATA.RDS.CForUL.CF,'sk','MarkerFaceColor','k','DisplayName','Lines 20260124')
xlabel("Frequency [Hz]")
ylabel("CF [m/T]")
title("CEB Magnetic Coupling")
adjustAxisLimits(0.05,'X')
adjustAxisLimits(0.05,'Y')
legend()
if saveplots
    exportgraphics(gcf, sprintf("MAG_CF_plots_Metatron/CEB_CF_Magnetic_Lines_20260124_plus_Sweeps.png"), 'Resolution', resolution);
end

%% Plot WEB Lines new + all Sweeps
set_plot_default_properties();
% close all
figure()
hold on
for i=1:numel(colorForPercentiles)
    plot(WEB_percentiles.freq, WEB_percentiles.percCF(i,:),'.','Color',colorForPercentiles(i),'DisplayName',sprintf("%d %% Sweep percentile",percToPlot(i)))
end
plot(WEB_CForUL_20251211.freq,WEB_CForUL_20251211.CF,'sk','MarkerFaceColor','k','DisplayName','Lines 20251211')
xlabel("Frequency [Hz]")
ylabel("CF [m/T]")
title("WEB Magnetic Coupling")
adjustAxisLimits(0.05,'X')
adjustAxisLimits(0.05,'Y')
legend()
if saveplots
    exportgraphics(gcf, sprintf("MAG_CF_plots_Metatron/WEB_CF_Magnetic_Lines_20251211_plus_Sweeps.png"), 'Resolution', resolution);
end




