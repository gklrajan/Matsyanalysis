clc; clearvars; close all;

%load('/Institut Curie/Lab/Projects/Scripts/ZebranalysisSystem/');

root = '/Institut Curie/Lab/Projects/Scripts/ZebranalysisSystem';

video = [root '/2018-02-16_freeswim.avi'];
csv   = [root, '/2018-02-16_freeswim.csv'];

tailBeats = csvread('2018-02-16_freeswim.csv');
tailBeats_corrected=[tailBeats(:,1),tailBeats(:,2)-tailBeats(1,2)];

tailBeats_corrected=[tailBeats(:,1),tailBeats(:,2)-tailBeats(1,2)];

% angle correction and filtering

tmp_delta_ori_filtered = tailBeats_corrected(:,2);
tmp_delta_ori_filtered(isnan(tmp_delta_ori_filtered))=0; % remoivng nans for filtering orientation
windowWidth = 31;
polynomialOrder = 3;
tmp_delta_ori_filtered = sgolayfilt(tmp_delta_ori_filtered, polynomialOrder, windowWidth);

%tailBeats = cell2mat(tailBeats);

% % VERGENCE
% 
% Langl                   = rad2deg(tmp_data(:,2));
% Langl(Langl>60)         = nan;
% Langl(Langl<-60)        = nan;
% 
% Rangl                   = rad2deg(tmp_data(:,3));
% Rangl(Rangl>60)         = nan;
% Rangl(Rangl<-60)        = nan;
% 
% tmp_verg                = Langl-Rangl;
% 
% tmp_Vfb          = fbtrim(tmp_verg,[0.5, 0.6, 0.01, 5, 5]); % Noise reduction of the differntiated vergence signal
% tmp_vbad         = isnan(tmp_Vfb);
% tmp_vlength      = length(tmp_Vfb);
% tmp_xraw         = 1:tmp_vlength;
% tmp_vcrop        = tmp_Vfb;
% tmp_vcrop(tmp_vbad)  = [];
% tmp_xraw(tmp_vbad)   = [];
% 
% tmp_Vint         = interp1(tmp_xraw, tmp_vcrop, 1:tmp_vlength);
% tmp_allbad       = isnan(tmp_Vint);
% tmp_Vint(tmp_allbad) = 0;
% 
% filter1 = ones(1,5)./5;   % low-pass filter with high cut-off frequency
% filter2 = ones(1,20)./20; % low-pass filter with low cut-off frequency
% 
% tmp_vergF1           = filtfilt(filter1, 1, tmp_Vint);
% tmp_vergF2           = filtfilt(filter2, 1, tmp_Vint);
% tmp_vergF1(tmp_vbad) = nan;
% tmp_vergF2(tmp_vbad) = nan;

% % ORIENTATION
% 
% tmp_ori              = rad2deg(tmp_data(:,1)); % fish orientation (left turn = decrease, right turn = increase)
% 
% % BOUTS
% 
% xpos = tmp_data(:,4); ypos = tmp_data(:,5);
% 
% dx = [0; diff(xpos)]; % distance between two consecutive x-coordinates
% dy = [0; diff(ypos)]; % distance between two consecutive y-coordinates
% dx(isnan(dx)) = 0;
% dy(isnan(dy)) = 0;
% 
% filterV = ones(1,10)./10;   % narrow filter for peak velocity/timing detection, captures better bout on
% 
% dxV = filtfilt(filterV, 1, dx);
% dyV = filtfilt(filterV, 1, dy);
% tmp_fdistV = sqrt(dxV.^2 + dyV.^2); % distance moved between iterations, in pixels
% tmp_fvelV  = tmp_fdistV.*100;  % convert to velocity in mm/s


% % cropping: interpolate the NaNs in the x- and y-position-vectors and filter them
% 
% sq_side = 100;
% 
% xbad        = isnan(xpos);  ybad        = isnan(ypos);
% xlength     = length(xpos); ylength     = length(ypos);
% xraw        = 1:xlength;    yraw        = 1:ylength;
% xcrop       = xpos;         ycrop       = ypos;
% xcrop(xbad) = [];           ycrop(ybad) = [];
% xraw(xbad)  = [];           yraw(ybad)  = [];
% 
% Xint          = (interp1(xraw, xcrop, 1:xlength)); Yint          = (interp1(yraw, ycrop, 1:ylength));
% allbadx       = isnan(Xint);                       allbady       = isnan(Yint);
% Xint(allbadx) = 0;                                 Yint(allbady) = 0;
% 
% filter = ones(1,5)./5;   % low-pass filter to supress jitter
% XintF  = round(filtfilt(filter, 1, Xint)); YintF  = round(filtfilt(filter, 1, Yint));
% 
% clear xbad ybad xraw yraw xcrop ycrop xlength ylength allbadx allbady x_ y_ Xint Yint


IM_          = [];

frame_start = 1;
frame_end   = 1200;

ww = frame_start:frame_end; ww = ww-ww(1);
% xx = tmp_vergF2(frame_start:frame_end);
% yy = tmp_ori(frame_start:frame_end);
zz = tmp_delta_ori_filtered(frame_start:frame_end);

obj    = mmread(video, frame_start:frame_end, [], false, true);

x_ = XintF(frame_start:frame_end); % restrict the x- and y-positions to the frames in each chunk of video
y_ = YintF(frame_start:frame_end);

IM = zeros(2*sq_side+1, 2*sq_side+1, size(obj.frames,2) ,'uint8'); % pre-allocation

%%
clc;
close all;

time            = ww/300; % time in seconds

axiscolor       = [1,1,1];
backgroundcolor = [1,1,0]; % color choice is a bit problematic if you want to set the background to transparent afterwards

set(0,'DefaultFigureColor', backgroundcolor);

figure('units','normalized','outerposition',[0 0 1 1]);

for ii=1:length(time)
    
    tmp_ = rgb2gray(obj.frames(ii).cdata);
    
    % copy the image-frame into larger frame of grey background to prevent violations of the boundary conditions
    tmp = 0.5*ones(obj.width + 2*sq_side, obj.height + 2*sq_side,'uint8');
    tmp(1+sq_side:obj.height+sq_side, 1+sq_side:obj.width+sq_side) = tmp_;
    
    try
        tmp = tmp(sq_side + y_(ii)-sq_side:sq_side+y_(ii)+sq_side, sq_side + x_(ii)-sq_side:sq_side+x_(ii)+sq_side);
    catch
        if x_(ii) == 0 || y_(ii) == 0
            % substitute black image-frame if the x- and y-positions
            % could not be determined
            tmp = zeros(2*sq_side+1, 2*sq_side+1 ,'uint8');
            
        end
    end
    
%     % vergences
%     ax1 = subplot(3,4,[1,2]);
%     plot(time(1:ii), xx(1:ii),'Color',[0,1,0.98], 'LineWidth', 3);
%     xlim([-0.3 max(time)]);
%     ylim([0 1.20*max(xx)]);
%     ax1.XColor     = 'none';
%     ax1.Color      = 'none';
%     ax1.Box        = 'off';
%     ax1.LineWidth  = 4;
%     ax1.FontSize   = 22; % its not possible to separately define size for ticklabels and axislables
%     ax1.FontWeight = 'Bold';
%     ax1.YColor     = axiscolor;
%     
%     text(-0., 1.20*max(xx),'Eye Vergence in deg','Color',[0,1,0.98],'FontSize',28,'FontWeight','Bold');
%     %axis off
    
    % swim velocity
    ax2 = subplot(3,4,[5,6]);
    plot(time(1:ii), zz(1:ii),'Color',[0,0.44, 0.73], 'LineWidth', 3);
    xlim([-0.3 max(time)]);
    ylim([0 1.05*max(zz)]);
    ax2.XColor     = 'none';
    ax2.Color      = 'none';
    ax2.Box        = 'off';
    ax2.LineWidth  = 4;
    ax2.FontSize   = 22; % its not possible to separately define size for ticklabels and axislables
    ax2.FontWeight = 'Bold';
    ax2.YColor     = axiscolor;
    
    text(-0,1.05*max(zz),'tail angle in deg','Color',[0,0.44, 0.73],'FontSize',28,'FontWeight','Bold');
    %axis off
    
%     % orientation
%     ax3 = subplot(3,4,[9,10]);
%     plot(time(1:ii), yy(1:ii),'Color',[1,0, 0], 'LineWidth', 3);
%     xlim([-0.3 max(time)]);
%     ylim([1.05*min(yy) 1.05*max(yy)]);
%     
%     ax3.Color = 'none';
%     ax3.XColor            = axiscolor;
%     ax3.LineWidth         = 4;
%     ax3.FontSize          = 22; % its not possible to separately define size for ticklabels and axislables
%     ax3.FontWeight        = 'Bold';
%     ax3.XLabel.String     = 'time in[s]';
%     ax3.XLabel.FontName   = 'Arial';
%     ax3.XLabel.FontSize   = 20;
%     ax3.XLabel.FontWeight = 'Bold';
%     ax3.YColor            = axiscolor;
%     ax3.Box               = 'off';
%     
%     text(-0., 1.05*max(yy),'Orientation in deg','Color',[1,0, 0],'FontSize',28,'FontWeight','Bold');
    
    % video
    ax4 = subplot(3,4,[3,4,7,8,11,12]);
    imagesc(tmp); colormap(gray);
    pbaspect(ax4,[1 1 1])
    axis off

    %pause(0.01);
    
    frame = getframe(gcf);
        
    if ii == 1
        [SIf, cm]    = rgb2ind(frame.cdata, 256);
        transp_color = double(SIf(1,1));
        
        imwrite(SIf, cm, [root,'/myfig.gif'], 'Loop',1, 'Delay',0, 'TransparentColor', transp_color); % check used colors before writing the image because of transparency!!!!
    elseif mod(ii,2)
        [SIf]        = rgb2ind(frame.cdata, cm);
        transp_color = double(SIf(1,1));
        
        imwrite(SIf, cm, [root,'/myfig.gif'], 'WriteMode','append', 'Delay',0, 'TransparentColor',  transp_color);
    end

end

%%
clear options
options.big     = true; % Use BigTIFF format
options.ask     = true;
options.message = true;
saveastiff(IM, [root '/cropped.tif'], options);

% or concatentate to IM_ matrix for later compression etc