%% Spheroid analysis
clear
warning off
clc
close all

level=2^8-1;
% Insert directory where the images are stored
datadir = 'C:\';
% Create a saving directory
savedir = [datadir, 'save_image\'];
if ~isfolder(savedir)
    mkdir(savedir)
end
% Read the files with initials 
files = dir([datadir, 'S1_*.tif']); %initials of the tiff files in sequence
%
cols=lines(numel(files)); % color to plot
% bg = double(imread([datadir, files(1).name]))/level;

clear spheroid

%%
N = numel(files);
time = NaN(1,N);
SpheroidArea = NaN(1, N);
H_mean = NaN(1, N);
spheroidCentroid = NaN(N, 2);
Ncell = NaN(1,N);
cellArea = cell(1,N);
cellCircularity = cell(1,N);
cellDistance = cell(1,N);
pWidth = [];
scale = 1.3; %edit scale according to your image resolution and pixel density. 
for k=1:numel(files)
    orig = double(imread([datadir, files(k).name]))/level;
    orig = rgb2gray(orig);
%     figure(1); clf
%     imshow(orig);
%     return
    I = imcomplement(orig); 
%     I1 = imsharpen(I,'Radius',2,'Amount',3);
%     I(I<0)=0;
    nbins=10;
    J = histeq(I);
    J1 = imadjustn(J,[0.8 0.95],[]);
%     figure(2); clf
%     imshow(J1);
    % Edge detection
    [~,threshold] = edge(J1,'sobel');
    fudgeFactor = 0.25; % play with this to get the edges properly
    BWs = edge(I,'sobel',threshold * fudgeFactor);
%     figure(3); clf
%     imshow(BWs);
    % Dilation of image
    se90 = strel('line',2,90);
    se0 = strel('line',2,0);

    BWsdil = imdilate(BWs,[se90 se0]);
%     figure(4); clf
%     imshow(BWsdil);
    BWdfill = imfill(BWsdil,'holes');
    BWerode = imerode(BWdfill, [se90 se0]);
%     figure(5); clf
%     imshow(BWerode);
    BW = BWerode; %imclearborder(BWerode);
    finalBW = bwareafilt(BW, [50 2e6], 4);
    finalBW = imfill(finalBW, 'holes');
    
    
    B = bwboundaries(finalBW);
    s = regionprops(finalBW, 'Area', 'Centroid', 'Perimeter', 'Circularity', 'MajorAxisLength', 'MinorAxisLength');
    
    Area = [s.Area];
    Circularity = [s.Circularity];
    
    % Store data
    index1 = find(Area==max(Area)); % spheroid index
    index2 = find(Area~=max(Area)); % cell indices

    SpheroidArea(k) = Area(index1); % spheroid area
    spheroidCentroid(k,:) = s(index1).Centroid;  % spheroid centroid
    
    cellCentroid = [s(index2).Centroid];
    cellCentroid = reshape(cellCentroid, 2,  numel(index2));
    cellCentroid = cellCentroid';

    cellArea{k} = Area(index2);
    cellCircularity{k} = Circularity(index2);
    SC = cellCentroid-spheroidCentroid(k,:);
    cellDistance{k} = sqrt(SC(:,1).^2+SC(:,2).^2);
    Ncell(k) =  numel(index2);
    
    fprintf('Number of cell found in frame %d is %d ...\n', k, numel(index2))
% return
     %% Cartesian to polar transform
    [ny,nx] = size(finalBW);
    imCentroid = [nx, ny]/2; % image centroids
    
    % displace the spheroid to center of image
    displaceSpheroid = imCentroid-spheroidCentroid(k,:); 
    shiftBW = imtranslate(finalBW, displaceSpheroid);
    % Remove single cells
    shiftBW = bwareafilt(shiftBW, [Area(index1)-1000 2e6], 4);
    % Transform the binary image into polar coordinates
    BW_pol = im_cart2pol(shiftBW); % binary image
    se = strel('disk', 5);  %EDIT
    BW_pol= imclose(BW_pol, se);  %EDIT
    % Find edges
    Ipoledge = edge(BW_pol);
    [nyp, nxp] = size(BW_pol);
    % Pixel to angular transform
    [row,col] = find(Ipoledge==1);
    xlin = linspace(0, 360,size(Ipoledge,2));
    ylin = linspace(0, size(Ipoledge,1), size(Ipoledge,1));
    ylin = ylin*scale;  %EDIT
    dispname = sprintf('t = %d s', k-1);
    time(k) = k-1;
    h= figure(1); clf
%     subplot(1,2,1)
    imshow(orig,[]); hold on
    cellfun(@(x) plot(x(:,2), x(:,1), '.' , 'Color', 'r','LineWidth',1.5,'DisplayName',dispname), B(index1))
    cellfun(@(x) plot(x(:,2), x(:,1), '.' , 'Color', 'g','LineWidth',1.5,'DisplayName',dispname), B(index2))
    plot(s(index1).Centroid(1), s(index1).Centroid(2), 'oy', 'MarkerFaceColor','y')
    if ~isempty(index2)
        plot(cellCentroid(:,1), cellCentroid(:,2), '.m', 'MarkerFaceColor','m')
    end
    text(50, 100, dispname, 'FontSize', 14, 'Color', 'y')
    saveas(gcf, [savedir, 'spheroid_', files(k).name])
    figure(2); clf
    imagesc(xlin, ylin,BW_pol); axis equal xy
    xlabel('$\theta$ [$\circ$]', 'Interpreter', 'latex');
    ylabel('$I$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 18)
    axis tight
    axis([0 360 0 nyp])
    
    colormap("gray")
    caxis([0 1])
    hold on
    plot(xlin(col), ylin(row), '.r', 'LineWidth',1.5)
    
 
    H_mean(k) = mean(ylin(row));
    H_std(k) = std(ylin(row));
    
    plot(xlin(col), H_mean(k)*ones(size(xlin(col))), '--y')

    xp = xlin(col);
    yp = ylin(row);
    ym = H_mean(k);
    [xp_new, yp_new] = dist_sort(xp, yp);
    plot(xp_new, yp_new, 'y')
    
    [W, H] = peak_analysis(xp_new,yp_new,ym);
    W_px = 2 *pi * ym * (W/360);
    pWidth{k} = W_px*scale; %in micron
    pHeight{k} = H*scale; %in micron
    
    saveas(gcf, [savedir, 'ploar_', files(k).name])
    save([savedir, 'peak_outputs.mat'], 'pWidth', 'pHeight') %this data is stored for peak height analysis
end

%%
figure(2); clf
subplot(1,3,1)
plot(time, sqrt(SpheroidArea/pi), 'om-', 'DisplayName', '$R_{s}$'); hold on
plot(time, H_mean, 'og-', 'DisplayName', '$H_{s}$'); 
xlabel('time [hr]', 'Interpreter','latex')
ylabel('Effective spheroid radius [px]', 'Interpreter','latex')
subplot(1,3,2)
plot(time, H_std, 'og-', 'DisplayName', '$H_{s}$'); 
xlabel('time [hr]', 'Interpreter','latex')
ylabel('Standard deviation[px]', 'Interpreter','latex')
L = legend;
set(L, 'Fontsize', 18, 'Interpreter', 'latex')
set(gca, 'FontSize', 18)
subplot(1,3,3)
plot(time, Ncell-Ncell(1), 'om-', 'DisplayName', '$R_{s}$'); hold on
xlabel('time [hr]', 'Interpreter','latex')
ylabel('Number of individual cells', 'Interpreter','latex')
set(gca, 'FontSize', 18)

%save directory file name
save([savedir, 'MV3_Ctrl_S1.mat'], 'SpheroidArea', 'spheroidCentroid', 'cellCentroid','cellArea','cellCircularity', 'cellDistance', 'Ncell', 'time', "xlin", 'ylin', 'H', 'col', 'row')


