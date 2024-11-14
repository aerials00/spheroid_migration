%% Spheroid analysis
clear
warning off
% clc
% close all

mainDirectory = 'C:\'; %insert main directory here 
ncases = 5; %adjust number based on number of spheroids 
nt = 25; %adjust time points for MV3 and A549
scale = 1.3; % um/px
Deff = NaN(ncases,nt);
for k=1:ncases
    subDirectory = [mainDirectory, sprintf('S%d', k), '\save_image\'];
    filename = dir([subDirectory, '*.mat']); %enter file name
    if ~isempty(filename)
        load([subDirectory, filename.name])
        Deff(k,:) = sqrt(SpheroidArea/pi)/sqrt(SpheroidArea(1)/pi);
%           Deff(k,:) = (SpheroidArea - SpheroidArea(1))/SpheroidArea(1); %sqrt(SpheroidArea/pi)/sqrt(SpheroidArea(1)/pi);
            
            ind10 = find(Deff(k,:)<1.1, 1, 'last');
             t_10(k) = time(ind10);
    end

% return
Deff_mean = nanmean(Deff,1);
Deff_std  = nanstd(Deff,1,1);
for k =1:ncases
    figure(2);
    plot(time, Deff(k,:)/Deff(k,1), 'o-', 'DisplayName', '$R_{s}$'); hold on
end
figure(2);
xlabel('Time [hr]', 'Interpreter','latex')
ylabel('Equivalent circular radius (-)', 'Interpreter','latex')
ylim ([0.9 2.4])
L = legend;
set(L, 'Fontsize', 20, 'Interpreter', 'latex')
set(gca, 'FontSize', 20)
errorbar(time, Deff_mean, Deff_std, 'ob-', 'MarkerFaceColor', 'b', 'LineWidth',1.5, 'DisplayName', '8mg/ml; MMP'); hold on
% errorbar(time, Deff_mean, Deff_std, 'diamondm-', 'MarkerFaceColor', 'w', 'LineWidth',1.5, 'DisplayName', '8 mg/ml; MMP'); hold on
save([mainDirectory, 'onset_time.mat'], 't_10'); %saves onset time for all spheroids 