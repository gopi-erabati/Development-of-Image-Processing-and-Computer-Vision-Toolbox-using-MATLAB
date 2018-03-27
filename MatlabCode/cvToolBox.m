function varargout = cvToolBox(varargin)
% CVTOOLBOX MATLAB code for cvToolBox.fig
%      CVTOOLBOX, by itself, creates a new CVTOOLBOX or raises the existing
%      singleton*.
%
%      H = CVTOOLBOX returns the handle to a new CVTOOLBOX or the handle to
%      the existing singleton*.
%
%      CVTOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CVTOOLBOX.M with the given input arguments.
%
%      CVTOOLBOX('Property','Value',...) creates a new CVTOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cvToolBox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cvToolBox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cvToolBox

% Last Modified by GUIDE v2.5 06-Jun-2017 19:05:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cvToolBox_OpeningFcn, ...
                   'gui_OutputFcn',  @cvToolBox_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cvToolBox is made visible.
function cvToolBox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cvToolBox (see VARARGIN)

% Choose default command line output for cvToolBox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%to set axes handles not to show ticks on axes
axes(handles.axes1);
set(gca,'XtickLabel',[],'YtickLabel',[]);

%to set axes handles not to show ticks on axes
axes(handles.axes2);
set(gca,'XtickLabel',[],'YtickLabel',[]);

set(handles.histogramButtonGroup,'visible','off');
set(handles.filterButtonGroup,'visible','off');
set(handles.edgesAndLinesButtonGroup, 'visible', 'off');
set(handles.contoursAndShapesButtonGroup, 'visible', 'off');
set(handles.featuresAndMatchesButtonGroup,'visible','off');
set(handles.cameraButtonGroup,'visible','off');

set(handles.histPanel,'visible','off');





% UIWAIT makes cvToolBox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cvToolBox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadImagePush.
function loadImagePush_Callback(hObject, eventdata, handles)
% hObject    handle to loadImagePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Function to get image
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.bmp;*.jpeg;*.png;*.gif','All Image Files';'*.*','All Files'}, 'Select an Image');
fileName = fullfile(pathname, filename);

% read image
img = imread(fileName);

axes(handles.axes1);
imshow(img);

handles.img = img;
guidata(hObject, handles);

set(handles.processPush,'Enable','on'); % set process button vissibility on


% --- Executes on button press in processPush.
function processPush_Callback(hObject, eventdata, handles)
% hObject    handle to processPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.saveText,'visible','off'); %set saved status label visible off

%no input images through general push button for calibration and stitching
if ((get(handles.calibrateRadio,'Value') == false) && (get(handles.stitchRadio,'Value') == false))
    % get the input image
    img = handles.img;
end

if ((get(handles.calibrateRadio,'Value') == true) || (get(handles.stitchRadio,'Value') == true))
    img = [];
end

%if check box is tick apply processing on result image!
if get(handles.resultCheck,'value') == true
    img = handles.img_out;
    axes(handles.axes1);
    imshow(img);
else
    if numel(size(img)) ~= 4
    axes(handles.axes1);
    imshow(img);
    end
end

if numel(size(img)) ~= 4  %to check for video or image

% for General tab
if get(handles.generalRadio, 'Value') == true
    
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    % For Noise
    if get(handles.noiseRadio,'Value') == true
        % for salt and pepper
        if get(handles.saltRadio, 'Value') == true
            %get the noise density value
            nDen = str2double(get(handles.saltDensityEdit, 'String'));
            img_out = imnoise(img,'salt & pepper', nDen); % add noise
            % show output image
            axes(handles.axes2);
            imshow(img_out);
            handles.img_out = img_out; guidata(hObject,handles);
        
        % for gaussian noise
        elseif get(handles.gaussianRadio, 'Value') == true
            %get mean and variance
            gMean = str2double(get(handles.meanEdit, 'String')); 
            gVar = str2double(get(handles.varianceEdit, 'String'));
            img_out = imnoise(img,'gaussian',gMean,gVar);
            % show output image
            axes(handles.axes2);
            imshow(img_out);
            handles.img_out = img_out; guidata(hObject,handles);
            
            %for speckl;e noise
        else
            sDen = str2double(get(handles.varianceSpeckleEdit, 'String')); 
            img_out = imnoise(img,'speckle',sDen);
                        % show output image
            axes(handles.axes2);
            imshow(img_out);
            handles.img_out = img_out; guidata(hObject,handles);
        end
        
        % for adding logo
    elseif get(handles.logoRadio, 'Value') == true
        %get data from lad push callback
        img_logo = handles.img_logo;
        
        if numel(size(img)) == 3
            
            % Get the size of the  logo - we need this
            % to properly mix the stuff in
            rows = size(img_logo, 1);
            cols = size(img_logo, 2);
            
            alpha = 0.4;
            alpha = repmat(alpha, rows, cols);
            %// Duplicate the alpha map to make this three channels - This is a colour image
            alphaResize = repmat(alpha, [1 1 3]);
            %// Make double to ensure mixing
            alphaResize = im2double(alphaResize);
            
            %make a duplicate of image
            img_out = img;
            
            % add logo at bottom right
            img_out(end-rows+1:end,end-cols+1:end,:) = alphaResize.*double(img_logo) + ...
                (1-alphaResize).*double(img_out(end-rows+1:end,end-cols+1:end,:));
            
            % show output image
            axes(handles.axes2);
            imshow(img_out);
            handles.img_out = img_out; guidata(hObject,handles);
            
        else
            if numel(size(img_logo)) == 3
                img_logo = rgb2gray(img_logo);
            end
            % Get the size of the  logo - we need this
            % to properly mix the stuff in
            rows = size(img_logo, 1);
            cols = size(img_logo, 2);
            
            alpha = 0.4;
            alpha = repmat(alpha, rows, cols);
            
            %// Make double to ensure mixing
            alphaResize = im2double(alpha);
            
            %make a duplicate of image
            img_out = img;
            
            % add logo at bottom right
            img_out(end-rows+1:end,end-cols+1:end) = alphaResize.*double(img_logo) + ...
                (1-alphaResize).*double(img_out(end-rows+1:end,end-cols+1:end));
            
            % show output image
            axes(handles.axes2);
            imshow(img_out);
            handles.img_out = img_out; guidata(hObject,handles);
            
        end
            
        %change color spaces
    elseif get(handles.colorspaceRadio,'Value') == true
        if get(handles.rgb2grayRadio, 'Value') == true
            img_out = rgb2gray(img);
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.rgb2hsvRadio, 'Value') == true
            img_out = rgb2hsv(img);
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.rgb2xyzRadio, 'Value') == true
            img_out = rgb2xyz(img);
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.rgb2ycbcrRadio, 'Value') == true
            img_out = rgb2ycbcr(img);
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.rgb2labRadio, 'Value') == true
            img_out = rgb2lab(img);
            handles.img_out = img_out; guidata(hObject,handles);
        else
            img_out = rgb2ntsc(img);
            handles.img_out = img_out; guidata(hObject,handles);
        end
        
        % show output image
        axes(handles.axes2);
        imshow(img_out);
        
    end
    
% for histogram TAB
elseif get(handles.histogramRadio, 'Value') == true
    % for histograms
    if get(handles.histogramRadio1, 'Value') == true
        set(handles.histPanel,'visible','on'); % to make histogram panel on
        % for plotting histograms
        if get(handles.plotHistogramRadio, 'Value') == true
            
            % for gray scale image
            if numel(size(img)) == 2 % gray image
                [counts,binLocations] = imhist(img); % create histogram
                ax1 = subplot(1,1,1,'Parent', handles.histPanel); % set axis and plot
                bar(ax1,binLocations,counts, 'BarWidth', 1);
                axis([0 255 0 max(counts)]);
            else % for rgb image
                [countsR,binLocationsR] = imhist(img(:,:,1)); %R
                [countsG,binLocationsG] = imhist(img(:,:,2));%G
                [countsB,binLocationsB] = imhist(img(:,:,3)); %B
                % show plots
                Ymax = max([countsR', countsG', countsB']);
                ax1 = subplot(3,1,1,'Parent', handles.histPanel);
                bar(ax1, binLocationsR, countsR, 'BarWidth', 1);
                axis([0 255 0 Ymax]);
                ax2 = subplot(3,1,2);
                bar(ax2, binLocationsG, countsG, 'BarWidth', 1);
                axis([0 255 0 Ymax]);
                ax3 = subplot(3,1,3);
                bar(ax3, binLocationsB, countsB, 'BarWidth', 1);
                axis([0 255 0 Ymax]);
            end
            
            %for equalisation of image
        else
            if numel(size(img)) == 3
                img = rgb2gray(img); %convert to gray
            end
            
            img_out = histeq(img); %equalise hist
            axes(handles.axes1)
            imshowpair(img, img_out, 'montage'); % show 2 images
            handles.img_out = img_out; guidata(hObject,handles);
            
            %show histogram
            [counts,binLocations] = imhist(img); % create histogram
            ax1 = subplot(2,1,1,'Parent', handles.histPanel); % set axis and plot
            bar(ax1,binLocations,counts, 'BarWidth', 1);
            axis([0 255 0 max(counts)]);
            
            [counts,binLocations] = imhist(img_out); % create histogram
            ax2 = subplot(2,1,2); % set axis and plot
            bar(ax2,binLocations,counts, 'BarWidth', 1);
            axis([0 255 0 max(counts)]);
            
            
        end
        
        % for morphology
    elseif get(handles.morphRadio,'Value') == true
        
        % for structuring element
        if get(handles.diamondRadio, 'Value') == true
            SE = strel('diamond', str2double(get(handles.radiusSEEdit,'String'))); % diamond SE
        elseif get(handles.diskRadio, 'Value') == true
            SE = strel('disk', str2double(get(handles.radiusSEEdit,'String'))); 
        else
            SE = strel('rectangle',[str2double(get(handles.lengthSEEdit,'String')), str2double(get(handles.widthSEEdit,'String'))]);
        end
        
        % checking for binary image
        if get(handles.binaryCheck, 'Value') == false
            
            %change image to gray scale
            if numel(size(img)) == 3
                img_temp = rgb2gray(img);
            else
                img_temp = img;
            end
            
            % for morph operatioons
            if get(handles.dilateRadio, 'Value') == true
                img_out = imdilate(img_temp,SE); % dilate image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.erodeRadio, 'Value') == true
                img_out = imerode(img_temp, SE); % erode image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.openRadio, 'Value') == true
                img_out = imopen(img_temp, SE); % open image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.closeRadio, 'Value') == true
                img_out = imclose(img_temp, SE); % open image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.hitmissRadio, 'Value') == true
                SE1 = strel('rectangle',[3 1]);
                SE2 = strel('disk', 1);
                img_out =  bwhitmiss(img_temp,SE1, SE2); % perform hit-miss
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.tophatRadio, 'Value') == true
                img_out = imtophat(img_temp, SE); %perofrm top hat filtering
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.bottomhatRadio, 'Value') == true
                img_out = imbothat(img_temp, SE); %perofrm bottom hat filtering
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.fillholesRadio, 'Value') == true
                img_out = imfill(img_temp,'holes');
                handles.img_out = img_out; guidata(hObject,handles);
            end
            
            % show images
            axes(handles.axes1)
            imshow(img_temp,[]);
            
            axes(handles.axes2)
            imshow(img_out, []);
            
        elseif get(handles.binaryCheck, 'value') == true
            
             if numel(size(img)) == 3
                img_temp = rgb2gray(img);
            else
                img_temp = img;
             end
            
            nTimes = str2double(get(handles.nTimesEdit, 'String'));
            
            % for morph operatioons
            if get(handles.dilateRadio, 'Value') == true
                img_out = imdilate(img_temp,ones(3)); % dilate image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.erodeRadio, 'Value') == true
                img_out = imerode(img_temp, ones(3)); % erode image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.openRadio, 'Value') == true
                img_out = bwmorph(img_temp, 'open', nTimes); % open image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.closeRadio, 'Value') == true
                img_out = bwmorph(img_temp, 'close', nTimes); % open image
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.hitmissRadio, 'Value') == true
                SE1 = strel('rectangle',[10 7]);
                SE2 = strel('disk', 6);
                img_out =  bwhitmiss(img_temp,SE1, SE2); % perform hit-miss
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.tophatRadio, 'Value') == true
                img_out = bwmorph(img_temp, 'tophat', nTimes); %perofrm top hat filtering
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.bottomhatRadio, 'Value') == true
                img_out = bwmorph(img_temp, 'bothat', nTimes); %perofrm bottom hat filtering
                handles.img_out = img_out; guidata(hObject,handles);
            elseif get(handles.fillholesRadio, 'Value') == true
                img_out = imfill(img_temp,'holes');
                handles.img_out = img_out; guidata(hObject,handles);
            end
            
            %show output
            axes(handles.axes2)
            imshow(img_out, []);
        end  
    end
    
% for filtering    
elseif get(handles.filterRadio, 'Value') == true
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    % for low pass filtering
    if get(handles.lpfRadio, 'Value') == true
        
        % for avearge filter
        if get(handles.averageRadio, 'Value') == true
            kSize = str2double(get(handles.kernelSizeEdit, 'String'));
            h = fspecial('average', kSize); % get avg  mask
            img_out = imfilter(img, h);
            handles.img_out = img_out; guidata(hObject,handles);
        
        % for gasusian filter
        elseif get(handles.gaussianFilterRadio, 'Value') == true
            gSigma = str2double(get(handles.gaussianSigmaEdit, 'String'));
            kSize = str2double(get(handles.kernelSizeEdit, 'String'));
            img_out = imgaussfilt(img, gSigma, 'FilterSize', kSize);
            handles.img_out = img_out; guidata(hObject,handles);
            
            %for median filter
        elseif get(handles.medianRadio, 'Value') == true
            if numel(size(img)) == 3
                img = rgb2gray(img);
                axes(handles.axes1); imshow(img,[]);
            end
            
            img_out = medfilt2(img); % median filter
            handles.img_out = img_out; guidata(hObject,handles);
          
            % for bilateral filter
        elseif get(handles.bilateralRadio, 'value') == true
            sSigma = str2double(get(handles.spatialSigmaEdit, 'String'));
            kSize = str2double(get(handles.kernelSizeEdit, 'String'));
            iSigma = str2double(get(handles.intensitySigmaEdit, 'String'));
            img_temp = (img./(max(max(max(img)))));
            img_out = bfilter2(double(img_temp), floor(kSize/2), [sSigma, iSigma]); % bilateral filter
            handles.img_out = img_out; guidata(hObject,handles);
        end
        axes(handles.axes2)
        imshow(img_out);
        
        % for high pass filtering
    elseif get(handles.hpfRadio, 'Value') == true
        
        % for sobel operator
        if get(handles.sobelRadio, 'value') == true
            h = fspecial('sobel'); 
            img_out = imfilter(img, h); % apply sobel operator
            handles.img_out = img_out; guidata(hObject,handles);
            
            %for prewitt operator
        elseif get(handles.prewittRadio, 'Value') == true
            h = fspecial('prewitt');
            img_out = imfilter(img, h); % apply sobel operator
            handles.img_out = img_out; guidata(hObject,handles);
            
            % for laplcaian of gaussian 
        elseif get(handles.logRadio, 'value') == true
            gSigma = str2double(get(handles.gaussianSigmaHpfEdit, 'String'));
            kSize = str2double(get(handles.kSizeHpfEdit, 'String'));            
            h = fspecial('log', kSize, gSigma);
            img_out = imfilter(img, h); % apply sobel operator
            handles.img_out = img_out; guidata(hObject,handles);
            
            % for image gradient
        elseif get(handles.gradientRadio, 'Value') == true
            
            % for gray scale
            img_temp = rgb2gray(img);
            [img_out,~] = imgradient(img_temp);
            handles.img_out = img_out; guidata(hObject,handles);
            
            %laplacian
        elseif get(handles.laplacianRadio, 'value') == true
            alphaLap = str2double(get(handles.alphaLapEdit, 'String'));
            h = fspecial('laplacian', alphaLap);
            img_out = imfilter(img, h); % apply laplacian 
            handles.img_out = img_out; guidata(hObject,handles);
            
            % for sharpen
        elseif get(handles.sharpenRadio, 'value') == true
            img_out = imsharpen(img); % to sharpen image
            handles.img_out = img_out; guidata(hObject,handles);
        end
        
        axes(handles.axes2)
        imshow(img_out);    
    end
    
    %for edges and lines
elseif get(handles.edgesAndLinesRadio, 'Value') == true
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    
    if numel(size(img)) == 3
    %convert to gray
        img_temp = rgb2gray(img);
    end
    
    %for edges
    if get(handles.edgeRadio,'Value') == true
        threshold = str2double(get(handles.thresholdEdgeEdit, 'String'));
        %for sobel edges
        if get(handles.sobelEdgeRadio, 'Value') == true
            img_out = edge(img_temp,'Sobel',threshold,'both'); % sobel edge detector
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.prewittEdgeRadio, 'Value') == true
            img_out = edge(img_temp,'Prewitt',threshold,'both'); % prewitt edge detector
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.robertsEdgeRadio, 'Value') == true
            img_out = edge(img_temp,'Roberts',threshold); % Roberts edge detector  
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.logEdgeRadio, 'Value') == true
            gSigma = str2double(get(handles.sigmaCannyEdit, 'String'));            
            img_out = edge(img_temp, 'log', threshold, gSigma); %LOG edges
            handles.img_out = img_out; guidata(hObject,handles);
        elseif get(handles.cannyEdgeRadio, 'Value') == true
            gSigma = str2double(get(handles.sigmaCannyEdit, 'String'));  
            thresholdHigh = str2double(get(handles.thresholdEdgeHighEdit, 'String'));
            img_out = edge(img_temp, 'Canny',[threshold thresholdHigh],gSigma);
            handles.img_out = img_out; guidata(hObject,handles);
            
        elseif get(handles.zeroCrossEdgeRadio, 'Value') == true
            h = fspecial('sobel');
            img_out = edge(img_temp, 'zerocross',threshold,h);
            handles.img_out = img_out; guidata(hObject,handles);
        end
        
        axes(handles.axes1)
        imshow(img_temp, []);
        
        axes(handles.axes2)
        imshow(img_out, []);  
        
        %for lines Hough transform
    elseif get(handles.linesRadio, 'Value') == true
        
        if numel(size(img)) == 3
        %convert to gray
            img_temp = rgb2gray(img);
        end
        
        % get the parameters
        threshold = str2double(get(handles.thresholdCannyLinesEdit, 'String'));
        thresholdH = str2double(get(handles.thresholdCannyLinesHighEdit, 'String'));
        minGap = str2double(get(handles.minGapLinesEdit, 'String'));
        minLength = str2double(get(handles.minLengthLinesEdit, 'String'));
        
        %detect edges using canny
        img_bw = edge(img_temp, 'Canny',[threshold thresholdH]);
               
        %apply hough transform
        [H,theta,rho] = hough(img_bw);
        
        % apply hough peaks
        P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
        
        % get hough lines
        lines = houghlines(img_bw,theta,rho,P,'FillGap',minGap,'MinLength',minLength);
        
        %show images
        axes(handles.axes1)
        imshow(img_temp, []);
        
        %draw lines on image
        axes(handles.axes2)
        imshow(img_temp), hold on
        max_len = 0;
        for k = 1:length(lines)
            xy = [lines(k).point1; lines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',1,'Color','green');
            
            % Plot beginnings and ends of lines
            plot(xy(1,1),xy(1,2),'x','LineWidth',0.5,'Color','yellow');
            plot(xy(2,1),xy(2,2),'x','LineWidth',0.5,'Color','red');
            
            % Determine the endpoints of the longest line segment
            len = norm(lines(k).point1 - lines(k).point2);
            if ( len > max_len)
                max_len = len;
                xy_long = xy;
            end
        end
        % highlight the longest line segment
        plot(xy_long(:,1),xy_long(:,2),'LineWidth',1,'Color','red');
        
        %for cirlces using hogh transform
    elseif get(handles.circlesRadio, 'Value') == true
        
        % Create waitbar.
        h = waitbar(0,'Applying hough to fing circles...');
        set(h,'Name','Hough Tranform progress');
        
        if numel(size(img)) == 3
            %convert to gray
            img_temp = rgb2gray(img);
        end
        
        %get parameters
        sensitivity = str2double(get(handles.sensitivityCirclesEdit, 'String'));
        threshold = str2double(get(handles.thresholdCirclesEdit, 'String'));
        minRadius = str2double(get(handles.minRadiusEdit, 'String'));
        maxRadius = str2double(get(handles.maxRadiusEdit, 'String'));
        strongestC = str2double(get(handles.strongestCirclesEdit, 'String'));
        waitbar(0.4);
        % get the hough circles
        [centers, radii, ~] = imfindcircles(img_temp,[minRadius maxRadius],'Sensitivity', sensitivity, 'EdgeThreshold', threshold);
        
        % check if strongest are avaialbel or not if not show all
        if size(centers,1) < strongestC
            strongestC = size(centers,1);
        end
        waitbar(0.8)
        %retain striongest circles according to metric values
        centersStrong = centers(1:strongestC,:);
        radiiStrong = radii(1:strongestC);
       
        %draw strongest circles
        axes(handles.axes1);
        imshow(img_temp, []);
        
        axes(handles.axes2);
        imshow(img_temp,[]);
        viscircles(centersStrong, radiiStrong,'EdgeColor','b');
        
        close(h); % close wait bar
        
    end        
    
    %for Contours and Shapes TAB
elseif get(handles.contoursAndShapesRadio, 'Value') == true
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    
    if numel(size(img)) == 3
        %convert to gray
        img_temp = rgb2gray(img);
    end
    
    %for contours
    if get(handles.drawContoursRadio, 'Value') == true

        % for contours
        [B,~] = bwboundaries(img_temp,'noholes');
        
        % show images
        axes(handles.axes1)
        imshow(img_temp, []);
        
        axes(handles.axes2);
        imshow(img_temp, []);
        hold on;
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
        end
        hold off;
        %for shapes
    elseif get(handles.drawShapesRadio, 'Value') == true
        
        % for bounding box
        if get(handles.boundingBoxRadio, 'Value') == true
            CC = bwconncomp(img_temp);  % get connected componenets
            st = regionprops(CC, 'BoundingBox' ); % create bounding box 
            
            axes(handles.axes1)
            imshow(img_temp, []);
        
            axes(handles.axes2)
            imshow(img_temp, []);
            hold on;
            for k = 1 : length(st)
                thisBB = st(k).BoundingBox;
                rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
                    'EdgeColor','r','LineWidth',2 )
            end
            
            hold off;
            
            % for circles bounding
        elseif get(handles.boundingCircleRadio, 'value') == true
            % for contours
            [B,~] = bwboundaries(img_temp,'noholes');
            
            % show images
            
            axes(handles.axes1)
            imshow(img_temp, []);
        
            axes(handles.axes2);
            imshow(img_temp, []);
            hold on;
            for k = 1:length(B)
                boundary = B{k};
                % for circle centre and radius
                [Xcenter, Ycenter, R] = SmallestEnclosingCircle(boundary(:,2)', boundary(:,1)');
                %draw circle
                viscircles([Xcenter, Ycenter], R,'EdgeColor','b');
            end
            hold off;
            
            % for convex hull
        elseif get(handles.convexHullRadio, 'value') == true

            img_out = bwconvhull(img_temp, 'objects');
            
            %show images
            axes(handles.axes1)
            imshow(img_temp, []);
        
            axes(handles.axes2)
            imshow(img_out);
            
            % for centroids 
        elseif get(handles.centroidRadio, 'Value') == true
            CC = bwconncomp(img_temp);  % get connected componenets
            
            %Calculate centroids for connected components in the image using regionprops.
            s = regionprops(CC,'centroid'); % region props for centroid
            
            % Concatenate structure array containing centroids into a single matrix.
            centroids = cat(1, s.Centroid);
            
            axes(handles.axes1)
            imshow(img_temp, []);
        
            % Display binary image with centroid locations superimposed.
            axes(handles.axes2)
            imshow(img_temp)
            hold on
            plot(centroids(:,1),centroids(:,2), 'b*')
            hold off
        end
        %for corners HARRIS
    elseif get(handles.cornersRadio, 'Value') == true
        sigmaC = str2double(get(handles.sigmaGaussianCornerEdit, 'String'));
        thresholdC = str2double(get(handles.thresholdCornerEdit, 'String'));
        radiusC = str2double(get(handles.radiusNMSEdit, 'String'));
        
        % from Peter toolbox
        [~, r, c] = harris(img_temp, sigmaC, thresholdC, radiusC, 0);
        
        axes(handles.axes1)
        imshow(img_temp, []);
        
        axes(handles.axes2)
        imshow(img_temp);
        hold on;
        plot(c,r,'yx');
        hold off;  
    end
    
    %for feature and Matches TAB
elseif get(handles.featuresAndMatchesRadio, 'Value') == true
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    
    if numel(size(img)) == 3
        %convert to gray
        img_temp = rgb2gray(img);
    end
    
    % for features extraction
    %for FAST
    if get(handles.fastRadio, 'value') == true
        %get parameters
        minQuality = str2double(get(handles.minQualityEdit, 'String'));
        minContrast = str2double(get(handles.minContrastEdit, 'String'));
        
        % detect FAST features
        points1 = detectFASTFeatures(img_temp,'MinQuality', minQuality,'MinContrast', minContrast);
             
        %show features in output
        %show input image
        axes(handles.axes1)
        imshow(img_temp, []);
        
        % show features on output image
        showFeatures = str2double(get(handles.showFeaturesEdit, 'String')); %get # of features to show
        
        axes(handles.axes2)
        imshow(img_temp); hold on;
        plot(points1.selectStrongest(showFeatures));
        
        %for Surf
    elseif get(handles.surfRadio, 'value') == true
        
        %get parameters
        thresholdSurf = str2double(get(handles.thresholdSurfEdit, 'String'));
        
        % detect SURF features
        points1 = detectSURFFeatures(img_temp,'MetricThreshold',thresholdSurf);
        
        %extract features
        [features1,validPoints1] = extractFeatures(img_temp, points1);
        
        %show features in output
        %show input image
        axes(handles.axes1)
        imshow(img_temp, []);
        
        % show features on output image
        showFeatures = str2double(get(handles.showFeaturesEdit, 'String')); %get # of features to show
        
%         figure('Name','Extracted features for Image 1');
        axes(handles.axes2)
        imshow(img_temp); hold on;
        plot(validPoints1.selectStrongest(showFeatures),'showOrientation',true);
        
        if get(handles.matchFeaturesRadio, 'Value') == true
            
            %get second image
            img_temp2 = handles.img2;
            
            
            % check whether you have second image
            
            if isempty(img_temp2)
                warndlg('No second Image','!! Warning !!');
            else
                if numel(size(img_temp2)) == 3
                    %convert to gray
                    img_temp2 = rgb2gray(img_temp2);
                end
                
                
                % detect FAST features for second image
                points2 = detectSURFFeatures(img_temp2,'MetricThreshold',thresholdSurf);
                
                %extract features
                [features2,validPoints2] = extractFeatures(img_temp2, points2);
                
                %show image1 with features
                axes(handles.axes1)
                imshow(img_temp,[]); hold on;
                plot(validPoints1.selectStrongest(showFeatures),'showOrientation',true);
                
                %show image2 with features
                axes(handles.axes2)
                imshow(img_temp2,[]); hold on;
                plot(validPoints2.selectStrongest(showFeatures),'showOrientation',true);
                               
                %get parameters
                matchThreshold = str2double(get(handles.matchThresholdEdit, 'String'));
                maxRatio = str2double(get(handles.maxRatioMatchEdit, 'String'));
                
                %match features
                indexPairs = matchFeatures(features1,features2, 'MatchThreshold', matchThreshold, 'MaxRatio', maxRatio);
                
                % Retrieve the locations of the corresponding points for each image.
                matchedPoints1 = validPoints1(indexPairs(:,1),:);
                matchedPoints2 = validPoints2(indexPairs(:,2),:);
                
                % Visualize the corresponding points
                %s hown in pop-up figure
                figure('Name','Matached Features');
                showMatchedFeatures(img_temp,img_temp2,matchedPoints1,matchedPoints2);
            end
        end
        
        %for HOG features
    elseif get(handles.hogRadio, 'Value') == true
        
        %extract hog features
        [~,hogVisualization] = extractHOGFeatures(img_temp);
        
        axes(handles.axes1)
        imshow(img_temp,[]);
        
        % show hog
        axes(handles.axes2);
        imshow(img_temp);
        hold on;
        plot(hogVisualization);
        
        % for SIFT
    elseif get(handles.siftRadio, 'Value') == true
        
         % Create waitbar.
        hsift = waitbar(0,'Applying SIFT to find matches...');
        set(hsift,'Name','SIFT in progress');
        
        % install vlfeat
        run('.\vlfeat\vlfeat-0.9.20\toolbox\vl_setup');
%         run('.\vlfeat\vlfeat-0.9.20\toolbox\vl_setup demo');
        
        img_temp = im2single(img_temp); % convert image to single
        
        %compute sift
        [f1,d1] = vl_sift(img_temp) ;
        
        showFeatures = str2double(get(handles.showFeaturesEdit, 'String')); %get # of features to show
        
        axes(handles.axes1)
        imshow(img_temp,[]);
        
        perm1 = randperm(size(f1,2)) ;
        sel1 = perm1(1:showFeatures) ;
        axes(handles.axes2)
        imshow(img_temp,[]); hold on;
        h11 = vl_plotframe(f1(:,sel1)) ;
        h21 = vl_plotframe(f1(:,sel1)) ;
        set(h11,'color','k','linewidth',3) ;
        set(h21,'color','y','linewidth',1) ;
        hold off;
        
        if get(handles.matchFeaturesRadio, 'Value') == true
            
            %get second image
            img_temp2 = handles.img2;
            
            
            % check whether you have second image
            
            if isempty(img_temp2)
                warndlg('No second Image','!! Warning !!');
            else
                if numel(size(img_temp2)) == 3
                    %convert to gray
                    img_temp2 = im2single(rgb2gray(img_temp2));
                end
                
                [f2, d2] = vl_sift(img_temp2) ;
                
                
                % show sift features for 2 images
                axes(handles.axes1)
                imshow(img_temp,[]); hold on;
                perm1 = randperm(size(f1,2)) ;
                sel1 = perm1(1:showFeatures) ;
                h11 = vl_plotframe(f1(:,sel1)) ;
                h21 = vl_plotframe(f1(:,sel1)) ;
                set(h11,'color','k','linewidth',3) ;
                set(h21,'color','y','linewidth',1) ;
                hold off;
                
                axes(handles.axes2)
                imshow(img_temp2,[]); hold on;
                perm2 = randperm(size(f2,2)) ;
                sel2 = perm2(1:showFeatures) ;
                h12 = vl_plotframe(f2(:,sel2)) ;
                h22 = vl_plotframe(f2(:,sel2)) ;
                set(h12,'color','k','linewidth',3) ;
                set(h22,'color','y','linewidth',1) ;
                hold off;
                
                
                [matches, scores] = vl_ubcmatch(d1, d2) ;
                
                % to show match lines
                [drop, perm] = sort(scores, 'descend') ;
                matches = matches(:, perm) ;
%                 scores  = scores(perm) ;
                
                matches = matches(:,1:showFeatures);
                
%                 figure('Name','Two images') ; clf ;
%                 imshow(cat(2, img_temp, img_temp2)) ;
%                 axis image off ;
                
                figure('Name','Matached Features') ; clf ;
                imshow(cat(2, img_temp, img_temp2)) ;
                
                xa = f1(1,matches(1,:)) ;
                xb = f2(1,matches(2,:)) + size(img_temp,2) ;
                ya = f1(2,matches(1,:)) ;
                yb = f2(2,matches(2,:)) ;
                
                hold on ;
                h = line([xa ; xb], [ya ; yb]) ;
                set(h,'linewidth', 1, 'color', 'b') ;
                
                vl_plotframe(f1(:,matches(1,:))) ;
                f2(1,:) = f2(1,:) + size(img_temp,2) ;
                vl_plotframe(f2(:,matches(2,:))) ;
                axis image off ;
                 
                 close(hsift)
            end
        end
    end
    
    %for CAMERA TAB
elseif get(handles.cameraRadio, 'Value') == true
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    
    %for calibration and undistortion
    if get(handles.calibrateRadio,'Value') == true
        
        h = waitbar(0,'calibrating in progress..');
        imgNames = handles.imgNames;
        
        %detect checkboard ppints
        [imagePoints,boardSize,~] = detectCheckerboardPoints(imgNames);
        
        % Generate world coordinates of the corners of the squares. Square size is in millimeters
        squareSize = 29;
        worldPoints = generateCheckerboardPoints(boardSize,squareSize);
        waitbar(0.4);
        % Calibrate the camera
        cameraParams = estimateCameraParameters(imagePoints,worldPoints);
        
         %detect checkboard ppints for one image
        [imagePoints,boardSize,~] = detectCheckerboardPoints(imgNames{10});
        
        % Generate world coordinates of the corners of the squares. Square size is in millimeters
        squareSize = 29;
        worldPoints = generateCheckerboardPoints(boardSize,squareSize);
        waitbar(0.8)
        %for extrinsic parameters
        [rotationMatrix,translationVector] = extrinsics(imagePoints,worldPoints,cameraParams);
        
        % for camera parameters
        camMatrix = cameraMatrix(cameraParams,rotationMatrix,translationVector);
        
        save('CameraParameters\camParameters.mat','camMatrix','rotationMatrix','translationVector');
        waitbar(0.9)
        % Remove lens distortion and display results.
        I = imread(imgNames{10});
        J1 = undistortImage(I,cameraParams);
        
        axes(handles.axes1)
        imshow(I);
        
        axes(handles.axes2)
        imshow(J1);
        
        figure; showExtrinsics(cameraParams);
        title('Extrinsics');
        waitbar(1)
        close(h)
        %for FUNDAMENTAL MATRIX using SIFT Matches
    elseif get(handles.fundamentalRadio,'Value') == true
        
        if numel(size(img)) == 3
            %convert to gray
            img_temp = rgb2gray(img);
        end
        
        axes(handles.axes1)
        imshow(img_temp,[]);
    
%          % Create waitbar.
%         hsift = waitbar(0,'Applying SIFT to find matches...');
%         set(hsift,'Name','SIFT in progress');
        
        % install vlfeat
        run('.\vlfeat\vlfeat-0.9.20\toolbox\vl_setup');
        
        img_temp = im2single(img_temp); % convert image to single
        
        %compute sift
        [f1,d1] = vl_sift(img_temp) ;
        
        %get second image
        img_temp2 = handles.img2;
        
        if numel(size(img_temp2)) == 3
            %convert to gray
            img_temp2 = im2single(rgb2gray(img_temp2));
        end
                %compute sift
        [f2, d2] = vl_sift(img_temp2) ;
        
        %find matches
        [matches, scores] = vl_ubcmatch(d1, d2) ;
        
        % to sort matches according to sxores
        [drop, perm] = sort(scores, 'descend') ;
        matches = matches(:, perm) ;
        
        %2d points of matched
        xa = f1(1,matches(1,:)) ;
        xb = f2(1,matches(2,:)) ;
        ya = f1(2,matches(1,:)) ;
        yb = f2(2,matches(2,:)) ;
        
        % Estimate the fundamental matrix from the two images
        
        points2d = [xa' ya' xb' yb'];
        
        [M,T1,T2]=normalHartley(points2d'); %hartley normalizatiomn of points
        
        F7p=funmat7p(M(:,1:7));
        F7p=T1'*F7p*T2;
        F7p=F7p./norm(F7p);
        
               
        %set visibility on to show F
        set(handles.fText,'visible','on');
        set(handles.f1Edit,'visible','on');
        set(handles.f2Edit,'visible','on');
        set(handles.f3Edit,'visible','on');
        set(handles.f4Edit,'visible','on');
        set(handles.f5Edit,'visible','on');
        set(handles.f6Edit,'visible','on');
        set(handles.f7Edit,'visible','on');
        set(handles.f8Edit,'visible','on');
        set(handles.f9Edit,'visible','on');
        
        set(handles.f1Edit,'String', num2str(F7p(1,1)));
        set(handles.f2Edit,'String', num2str(F7p(1,2)));
        set(handles.f3Edit,'String', num2str(F7p(1,3)));
        set(handles.f4Edit,'String', num2str(F7p(2,1)));
        set(handles.f5Edit,'String', num2str(F7p(2,2)));
        set(handles.f6Edit,'String', num2str(F7p(2,3)));
        set(handles.f7Edit,'String', num2str(F7p(3,1)));
        set(handles.f8Edit,'String', num2str(F7p(3,2)));
        set(handles.f9Edit,'String', num2str(F7p(3,3)));
        
        save('CameraParameters\fundamentalMatrix.mat','F7p');
        
        %for epipolar lines
    elseif get(handles.epipolarRadio,'Value') == true
        
        if numel(size(img)) == 3
            %convert to gray
            img_temp = rgb2gray(img);
        end
        h = waitbar(0,'Drawing epipolar lines...pls wait');
        
        %matching and finding fundamental matrix using SURF
        % detect SURF features
        points1 = detectSURFFeatures(img_temp);
        
        %extract features
        [features1,validPoints1] = extractFeatures(img_temp, points1);
%         validPoints1 = validPoints1.selectStrongest(20);

        
        %get second image
        img_temp2 = handles.img2;    
        
        if numel(size(img_temp2)) == 3
            %convert to gray
            img_temp2 = rgb2gray(img_temp2);
        end
        
        
        % detect FAST features for second image
        points2 = detectSURFFeatures(img_temp2);
        
        %extract features
        [features2,validPoints2] = extractFeatures(img_temp2, points2);
%         validPoints2 = validPoints2.selectStrongest(20);

        waitbar(0.5);

        %match features
        indexPairs = matchFeatures(features1,features2);
        
        % Retrieve the locations of the corresponding points for each image.
        matchedPoints1 = validPoints1(indexPairs(:,1),:);
        matchedPoints2 = validPoints2(indexPairs(:,2),:);
        
        [fLMedS,inliers] = estimateFundamentalMatrix(matchedPoints1,...
            matchedPoints2,'NumTrials',4000);
        
        % CONVERSIONS REQUIRED
        %%%% conversion of matched_points1 to matrix
        sS = 20;
        xm1 = matchedPoints1.selectStrongest(sS);
        m1 = zeros(sS,2,'double');
        newInliners = logical(sS);
        wielCP = size(matchedPoints1);
        wielCP = wielCP(1,1);
        
        for y = 1:sS
            am1 = xm1(y);
            m1(y,1) = am1.Location(1,1);
            m1(y,2) = am1.Location(1,2);
            for z = 1 : wielCP
                b = matchedPoints1(z);
                temp_a1 = am1.Location(1,1);
                temp_a2 = am1.Location(1,2);
                temp_b1 = b.Location(1,1);
                temp_b2 = b.Location(1,2);
                if temp_a1 == temp_b1 && temp_a2 == temp_b2
                    newInliners(y) = inliers(z);
                end
            end
        end
        newInliners = newInliners';
        %%%% conversion of matched_points2 to matrix
        xm2 = matchedPoints2.selectStrongest(sS);
        m2 = zeros(sS,2,'double');
        for y = 1:sS
            am2 = xm2(y);
            m2(y,1) = am2.Location(1,1);
            m2(y,2) = am2.Location(1,2);
        end
        %%%%%%%%
        waitbar(0.75);
        
                %show features in output
        %show input image
        axes(handles.axes1)
        imshow(img_temp, []); hold on;
        plot(m1(newInliners,1), m1(newInliners,2), 'go');
        
        %draw epipolar lines in first image
        epiLines = epipolarLine(fLMedS', m2(newInliners, :));
        
%         Compute the intersection points of the lines and the image border.
        points = lineToBorderPoints(epiLines,size(img_temp));
        
        %         Show the epipolar lines in the first image
        line(points(:,[1,3])',points(:,[2,4])','LineWidth',1)
        hold off;
        
        
        %show image2 with features
        axes(handles.axes2)
        imshow(img_temp2,[]); hold on;
        plot(m2(newInliners,1), m2(newInliners,2), 'go');
        epiLines = epipolarLine(fLMedS, m1(newInliners, :)); %draw epipolar lines in 2 image
        points = lineToBorderPoints(epiLines,size(img_temp2)); %         Compute the intersection points of the lines and the image border.
        line(points(:,[1,3])',points(:,[2,4])', 'LineWidth',1) %         Show the epipolar lines in the 2 image
        waitbar(1);
        close(h);
        
        %for Homography and Mosaic
    elseif get(handles.homographyRadio,'Value') == true
        
     
        
        
        %get second image
        img_temp2 = handles.img2;
        

        h = waitbar(0,'Computing homography...');
        
        %calculate homography
        [ H,~ ] = homography( img_temp2, img);
        
        waitbar(0.5);
        
        %set visibility on to show F
        set(handles.hText,'visible','on');
        set(handles.f1Edit,'visible','on');
        set(handles.f2Edit,'visible','on');
        set(handles.f3Edit,'visible','on');
        set(handles.f4Edit,'visible','on');
        set(handles.f5Edit,'visible','on');
        set(handles.f6Edit,'visible','on');
        set(handles.f7Edit,'visible','on');
        set(handles.f8Edit,'visible','on');
        set(handles.f9Edit,'visible','on');
        
        set(handles.f1Edit,'String', num2str(H(1,1)));
        set(handles.f2Edit,'String', num2str(H(1,2)));
        set(handles.f3Edit,'String', num2str(H(1,3)));
        set(handles.f4Edit,'String', num2str(H(2,1)));
        set(handles.f5Edit,'String', num2str(H(2,2)));
        set(handles.f6Edit,'String', num2str(H(2,3)));
        set(handles.f7Edit,'String', num2str(H(3,1)));
        set(handles.f8Edit,'String', num2str(H(3,2)));
        set(handles.f9Edit,'String', num2str(H(3,3)));
        
        save('CameraParameters\homographyMatrix.mat','H');
        
        waitbar(0.75);
        
        tform = maketform('projective',H');
        img21 = imtransform(img_temp2,tform); % reproject img2
        
        [M1 N1 dim] = size(img);
        [M2 N2 ~] = size(img_temp2);
        
        % do the mosaic
        pt = zeros(3,4);
        pt(:,1) = H*[1;1;1];
        pt(:,2) = H*[N2;1;1];
        pt(:,3) = H*[N2;M2;1];
        pt(:,4) = H*[1;M2;1];
        x2 = pt(1,:)./pt(3,:);
        y2 = pt(2,:)./pt(3,:);
        
        up = round(min(y2));
        Yoffset = 0;
        if up <= 0
            Yoffset = -up+1;
            up = 1;
        end
        
        left = round(min(x2));
        Xoffset = 0;
        if left<=0
            Xoffset = -left+1;
            left = 1;
        end
        
        [M3,N3,~] = size(img21);
        imgout(up:up+M3-1,left:left+N3-1,:) = img21;
        % img1 is above img21
        imgout(Yoffset+1:Yoffset+M1,Xoffset+1:Xoffset+N1,:) = img;
        waitbar(1);
        close(h)
        figure, imshow(imgout);
        
        %for stictching
    elseif get(handles.stitchRadio,'Value') == true
        
        h = waitbar(0,'Stitching...');
        
        imgS = handles.imgS;
        
        for i = 1:length(imgS)-1
       
        axes(handles.axes1); imshow(imgS{i});
        axes(handles.axes2); imshow(imgS{i+1});
        imgout = imStitch( imgS{i}, imgS{i+1});
        imgS{i+1} = imgout;
        waitbar(1/length(imgS))
        end
        
        waitbar(1)
        close(h)
        figure, imshow(imgout);
        
        
    end
    
end

else
    
    set(handles.histPanel,'visible','off'); % to make histogram panel off
    if get(handles.generalRadio, 'Value') == true
    
    
    % For Noise
    if get(handles.noiseRadio,'Value') == true
        % for salt and pepper
        if get(handles.saltRadio, 'Value') == true
            %get the noise density value
            nDen = str2double(get(handles.saltDensityEdit, 'String'));
           
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imnoise(img(:,:,:,i),'salt & pepper', nDen); % add noise
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
            
        
        % for gaussian noise
        elseif get(handles.gaussianRadio, 'Value') == true
            %get mean and variance
            gMean = str2double(get(handles.meanEdit, 'String')); 
            gVar = str2double(get(handles.varianceEdit, 'String'));
            
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imnoise(img(:,:,:,i),'gaussian',gMean,gVar);
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
            
            %for speckl;e noise
        else
            sDen = str2double(get(handles.varianceSpeckleEdit, 'String')); 
            
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imnoise(img(:,:,:,i),'speckle',sDen);
                        % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        end
        
               
        %change color spaces
    elseif get(handles.colorspaceRadio,'Value') == true
        if get(handles.rgb2grayRadio, 'Value') == true
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,i) = rgb2gray(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
            
        elseif get(handles.rgb2hsvRadio, 'Value') == true
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = rgb2hsv(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        elseif get(handles.rgb2xyzRadio, 'Value') == true
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = rgb2xyz(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        elseif get(handles.rgb2ycbcrRadio, 'Value') == true
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = rgb2ycbcr(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        elseif get(handles.rgb2labRadio, 'Value') == true
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = rgb2lab(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        else
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = rgb2ntsc(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        end
        
              
    end
    elseif get(handles.filterRadio, 'Value') == true
    
    % for low pass filtering
    if get(handles.lpfRadio, 'Value') == true
        
        % for avearge filter
        if get(handles.averageRadio, 'Value') == true
            kSize = str2double(get(handles.kernelSizeEdit, 'String'));
            h = fspecial('average', kSize); % get avg  mask
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imfilter(img(:,:,:,i),h);
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                  
        % for gasusian filter
        elseif get(handles.gaussianFilterRadio, 'Value') == true
            gSigma = str2double(get(handles.gaussianSigmaEdit, 'String'));
            kSize = str2double(get(handles.kernelSizeEdit, 'String'));
            
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imgaussfilt(img(:,:,:,i), gSigma, 'FilterSize', kSize);
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                        
            %for median filter
        elseif get(handles.medianRadio, 'Value') == true
                 
             for i = 1:size(img,4)
                imshow(rgb2gray(img(:,:,:,i)),'Parent',handles.axes1);
                img_out(:,:,i) = medfilt2(rgb2gray(img(:,:,:,i)));
                % show output image
                imshow(img_out(:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
            
                 
            % for bilateral filter
        elseif get(handles.bilateralRadio, 'value') == true
            sSigma = str2double(get(handles.spatialSigmaEdit, 'String'));
            kSize = str2double(get(handles.kernelSizeEdit, 'String'));
            iSigma = str2double(get(handles.intensitySigmaEdit, 'String'));
            
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_temp = (img(:,:,:,i)./(max(max(max(img(:,:,:,i))))));
                img_out(:,:,:,i) = bfilter2(double(img_temp), floor(kSize/2), [sSigma, iSigma]); % bilateral filter
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
        end
        
        % for high pass filtering
    elseif get(handles.hpfRadio, 'Value') == true
        
        % for sobel operator
        if get(handles.sobelRadio, 'value') == true
            h = fspecial('sobel'); 
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imfilter(img(:,:,:,i),h);
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                        
            %for prewitt operator
        elseif get(handles.prewittRadio, 'Value') == true
            h = fspecial('prewitt');
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imfilter(img(:,:,:,i),h);
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                       
            % for laplcaian of gaussian 
        elseif get(handles.logRadio, 'value') == true
            gSigma = str2double(get(handles.gaussianSigmaHpfEdit, 'String'));
            kSize = str2double(get(handles.kSizeHpfEdit, 'String'));            
            h = fspecial('log', kSize, gSigma);
            for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imfilter(img(:,:,:,i),h);
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                        
            % for image gradient
        elseif get(handles.gradientRadio, 'Value') == true
            
            for i = 1:size(img,4)
                imshow(rgb2gray(img(:,:,:,i)),'Parent',handles.axes1);
                img_temp = rgb2gray(img(:,:,:,i));
                img_out(:,:,i) = imgradient(img_temp);
                % show output image
                imshow(img_out(:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                       
            %laplacian
        elseif get(handles.laplacianRadio, 'value') == true
            alphaLap = str2double(get(handles.alphaLapEdit, 'String'));
            h = fspecial('laplacian', alphaLap);
             for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imfilter(img(:,:,:,i),h);
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
            end
                        
            % for sharpen
        elseif get(handles.sharpenRadio, 'value') == true
             for i = 1:size(img,4)
                imshow(img(:,:,:,i),'Parent',handles.axes1);
                img_out(:,:,:,i) = imsharpen(img(:,:,:,i));
                % show output image
                imshow(img_out(:,:,:,i),'Parent',handles.axes2); pause(0.000001);
             end
        end
        
          
    end
    end
end



    


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in noiseRadio.
function noiseRadio_Callback(hObject, eventdata, handles)
% hObject    handle to noiseRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noiseRadio


% --- Executes on button press in saltRadio.
function saltRadio_Callback(hObject, eventdata, handles)
% hObject    handle to saltRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saltRadio

% --- Executes on button press in gaussainRadio.
function gaussianRadio_Callback(hObject, eventdata, handles)
% hObject    handle to saltRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gaussianRadio

function saltDensityEdit_Callback(hObject, eventdata, handles)
% hObject    handle to saltDensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saltDensityEdit as text
%        str2double(get(hObject,'String')) returns contents of saltDensityEdit as a double


% --- Executes during object creation, after setting all properties.
function saltDensityEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saltDensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gaussianFilterRadio.
function gaussianFilterRadio_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianFilterRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gaussianFilterRadio



function meanEdit_Callback(hObject, eventdata, handles)
% hObject    handle to meanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meanEdit as text
%        str2double(get(hObject,'String')) returns contents of meanEdit as a double


% --- Executes during object creation, after setting all properties.
function meanEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function varianceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to varianceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varianceEdit as text
%        str2double(get(hObject,'String')) returns contents of varianceEdit as a double


% --- Executes during object creation, after setting all properties.
function varianceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varianceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in speckleRadio.
function speckleRadio_Callback(hObject, eventdata, handles)
% hObject    handle to speckleRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speckleRadio



function varianceSpeckleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to varianceSpeckleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varianceSpeckleEdit as text
%        str2double(get(hObject,'String')) returns contents of varianceSpeckleEdit as a double


% --- Executes during object creation, after setting all properties.
function varianceSpeckleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varianceSpeckleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in logoPush.
function logoPush_Callback(hObject, eventdata, handles)
% hObject    handle to logoPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Function to get image
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'}, 'Select a Logo');
fileName = fullfile(pathname, filename);

% read image
[logo] = imread(fileName);

handles.img_logo = logo;
guidata(hObject, handles);


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag')
    case 'generalRadio'
        set(handles.histogramButtonGroup,'visible','off');
        set(handles.generalButtonGroup,'visible','on');
        set(handles.filterButtonGroup,'visible','off');
        set(handles.edgesAndLinesButtonGroup,'visible','off');
        set(handles.contoursAndShapesButtonGroup,'visible','off');
        set(handles.featuresAndMatchesButtonGroup,'visible','off');
        set(handles.cameraButtonGroup,'visible','off');
    case 'histogramRadio'
        set(handles.generalButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'visible','on');
        set(handles.filterButtonGroup,'visible','off');
         set(handles.edgesAndLinesButtonGroup,'visible','off');
         set(handles.contoursAndShapesButtonGroup,'visible','off');
         set(handles.featuresAndMatchesButtonGroup,'visible','off');
        set(handles.cameraButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'position',get(handles.generalButtonGroup,'position'));
    case 'filterRadio'
        set(handles.generalButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'visible','off');
        set(handles.filterButtonGroup,'visible','on');
         set(handles.edgesAndLinesButtonGroup,'visible','off');
         set(handles.contoursAndShapesButtonGroup,'visible','off');
         set(handles.featuresAndMatchesButtonGroup,'visible','off');
         set(handles.cameraButtonGroup,'visible','off');
        set(handles.filterButtonGroup,'position',get(handles.generalButtonGroup,'position'));
    case 'edgesAndLinesRadio'
        set(handles.generalButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'visible','off');
        set(handles.filterButtonGroup,'visible','off');
        set(handles.edgesAndLinesButtonGroup,'visible','on');  
        set(handles.contoursAndShapesButtonGroup,'visible','off');
        set(handles.featuresAndMatchesButtonGroup,'visible','off');
        set(handles.cameraButtonGroup,'visible','off');
        set(handles.edgesAndLinesButtonGroup,'position',get(handles.generalButtonGroup,'position'));
    case 'contoursAndShapesRadio'
        set(handles.generalButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'visible','off');
        set(handles.filterButtonGroup,'visible','off');
        set(handles.edgesAndLinesButtonGroup,'visible','off');       
        set(handles.contoursAndShapesButtonGroup,'visible','on');
        set(handles.featuresAndMatchesButtonGroup,'visible','off');
        set(handles.cameraButtonGroup,'visible','off');
        set(handles.contoursAndShapesButtonGroup,'position',get(handles.generalButtonGroup,'position'));
    case 'featuresAndMatchesRadio'
        set(handles.generalButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'visible','off');
        set(handles.filterButtonGroup,'visible','off');
        set(handles.edgesAndLinesButtonGroup,'visible','off');       
        set(handles.contoursAndShapesButtonGroup,'visible','off');
        set(handles.featuresAndMatchesButtonGroup,'visible','on');
        set(handles.cameraButtonGroup,'visible','off');
        set(handles.featuresAndMatchesButtonGroup,'position',get(handles.generalButtonGroup,'position'));
    case 'cameraRadio'
        set(handles.generalButtonGroup,'visible','off');
        set(handles.histogramButtonGroup,'visible','off');
        set(handles.filterButtonGroup,'visible','off');
        set(handles.edgesAndLinesButtonGroup,'visible','off');       
        set(handles.contoursAndShapesButtonGroup,'visible','off');
        set(handles.featuresAndMatchesButtonGroup,'visible','off');
        set(handles.cameraButtonGroup,'visible','on');
        set(handles.cameraButtonGroup,'position',get(handles.generalButtonGroup,'position'));
end
         


% --- Executes when selected object is changed in histogramButtonGroup.
function histogramButtonGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in histogramButtonGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'histogramRadio1'
            set(handles.histPanel,'visible','on');
    case 'morphRadio'
        set(handles.histPanel,'visible','off');
end
        



function radiusSEEdit_Callback(hObject, eventdata, handles)
% hObject    handle to radiusSEEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radiusSEEdit as text
%        str2double(get(hObject,'String')) returns contents of radiusSEEdit as a double


% --- Executes during object creation, after setting all properties.
function radiusSEEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiusSEEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lengthSEEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lengthSEEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lengthSEEdit as text
%        str2double(get(hObject,'String')) returns contents of lengthSEEdit as a double


% --- Executes during object creation, after setting all properties.
function lengthSEEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lengthSEEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function widthSEEdit_Callback(hObject, eventdata, handles)
% hObject    handle to widthSEEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthSEEdit as text
%        str2double(get(hObject,'String')) returns contents of widthSEEdit as a double


% --- Executes during object creation, after setting all properties.
function widthSEEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to widthSEEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in binaryCheck.
function binaryCheck_Callback(hObject, eventdata, handles)
% hObject    handle to binaryCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of binaryCheck



function nTimesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nTimesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nTimesEdit as text
%        str2double(get(hObject,'String')) returns contents of nTimesEdit as a double


% --- Executes during object creation, after setting all properties.
function nTimesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nTimesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kernelSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to kernelSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernelSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of kernelSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function kernelSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernelSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gaussianSigmaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianSigmaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gaussianSigmaEdit as text
%        str2double(get(hObject,'String')) returns contents of gaussianSigmaEdit as a double


% --- Executes during object creation, after setting all properties.
function gaussianSigmaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gaussianSigmaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spatialSigmaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to spatialSigmaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spatialSigmaEdit as text
%        str2double(get(hObject,'String')) returns contents of spatialSigmaEdit as a double


% --- Executes during object creation, after setting all properties.
function spatialSigmaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatialSigmaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function intensitySigmaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to intensitySigmaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of intensitySigmaEdit as text
%        str2double(get(hObject,'String')) returns contents of intensitySigmaEdit as a double


% --- Executes during object creation, after setting all properties.
function intensitySigmaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intensitySigmaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kSizeHpfEdit_Callback(hObject, eventdata, handles)
% hObject    handle to kSizeHpfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kSizeHpfEdit as text
%        str2double(get(hObject,'String')) returns contents of kSizeHpfEdit as a double


% --- Executes during object creation, after setting all properties.
function kSizeHpfEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kSizeHpfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gaussianSigmaHpfEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianSigmaHpfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gaussianSigmaHpfEdit as text
%        str2double(get(hObject,'String')) returns contents of gaussianSigmaHpfEdit as a double


% --- Executes during object creation, after setting all properties.
function gaussianSigmaHpfEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gaussianSigmaHpfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphaLapEdit_Callback(hObject, eventdata, handles)
% hObject    handle to alphaLapEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaLapEdit as text
%        str2double(get(hObject,'String')) returns contents of alphaLapEdit as a double


% --- Executes during object creation, after setting all properties.
function alphaLapEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaLapEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdEdgeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdEdgeEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdEdgeEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdEdgeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaCannyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaCannyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaCannyEdit as text
%        str2double(get(hObject,'String')) returns contents of sigmaCannyEdit as a double


% --- Executes during object creation, after setting all properties.
function sigmaCannyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaCannyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdEdgeHighEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdEdgeHighEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdEdgeHighEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdEdgeHighEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdCannyLinesHighEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCannyLinesHighEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCannyLinesHighEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdCannyLinesHighEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdCannyLinesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCannyLinesEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCannyLinesEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdCannyLinesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minGapLinesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minGapLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minGapLinesEdit as text
%        str2double(get(hObject,'String')) returns contents of minGapLinesEdit as a double


% --- Executes during object creation, after setting all properties.
function minGapLinesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minGapLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minLengthLinesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minLengthLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minLengthLinesEdit as text
%        str2double(get(hObject,'String')) returns contents of minLengthLinesEdit as a double


% --- Executes during object creation, after setting all properties.
function minLengthLinesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLengthLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxRadiusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of maxRadiusEdit as a double


% --- Executes during object creation, after setting all properties.
function maxRadiusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minRadiusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of minRadiusEdit as a double


% --- Executes during object creation, after setting all properties.
function minRadiusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sensitivityCirclesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to sensitivityCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensitivityCirclesEdit as text
%        str2double(get(hObject,'String')) returns contents of sensitivityCirclesEdit as a double


% --- Executes during object creation, after setting all properties.
function sensitivityCirclesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensitivityCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdCirclesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCirclesEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCirclesEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdCirclesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function strongestCirclesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to strongestCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strongestCirclesEdit as text
%        str2double(get(hObject,'String')) returns contents of strongestCirclesEdit as a double


% --- Executes during object creation, after setting all properties.
function strongestCirclesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strongestCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdEdgeEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdEdgeEdit as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaCannyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaCannyEdit as text
%        str2double(get(hObject,'String')) returns contents of sigmaCannyEdit as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaCannyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdEdgeHighEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdEdgeHighEdit as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdEdgeHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCannyLinesHighEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCannyLinesHighEdit as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesHighEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCannyLinesEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCannyLinesEdit as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCannyLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to minGapLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minGapLinesEdit as text
%        str2double(get(hObject,'String')) returns contents of minGapLinesEdit as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minGapLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to minLengthLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minLengthLinesEdit as text
%        str2double(get(hObject,'String')) returns contents of minLengthLinesEdit as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLengthLinesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of maxRadiusEdit as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minRadiusEdit as text
%        str2double(get(hObject,'String')) returns contents of minRadiusEdit as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minRadiusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to sensitivityCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensitivityCirclesEdit as text
%        str2double(get(hObject,'String')) returns contents of sensitivityCirclesEdit as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensitivityCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCirclesEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCirclesEdit as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to strongestCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strongestCirclesEdit as text
%        str2double(get(hObject,'String')) returns contents of strongestCirclesEdit as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strongestCirclesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radiusNMSEdit_Callback(hObject, eventdata, handles)
% hObject    handle to radiusNMSEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radiusNMSEdit as text
%        str2double(get(hObject,'String')) returns contents of radiusNMSEdit as a double


% --- Executes during object creation, after setting all properties.
function radiusNMSEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiusNMSEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaGaussianCornerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaGaussianCornerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaGaussianCornerEdit as text
%        str2double(get(hObject,'String')) returns contents of sigmaGaussianCornerEdit as a double


% --- Executes during object creation, after setting all properties.
function sigmaGaussianCornerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaGaussianCornerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdCornerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCornerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdCornerEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdCornerEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdCornerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdCornerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minQualityEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minQualityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minQualityEdit as text
%        str2double(get(hObject,'String')) returns contents of minQualityEdit as a double


% --- Executes during object creation, after setting all properties.
function minQualityEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minQualityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minContrastEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minContrastEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minContrastEdit as text
%        str2double(get(hObject,'String')) returns contents of minContrastEdit as a double


% --- Executes during object creation, after setting all properties.
function minContrastEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minContrastEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresholdSurfEdit_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdSurfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdSurfEdit as text
%        str2double(get(hObject,'String')) returns contents of thresholdSurfEdit as a double


% --- Executes during object creation, after setting all properties.
function thresholdSurfEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdSurfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function showFeaturesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to showFeaturesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of showFeaturesEdit as text
%        str2double(get(hObject,'String')) returns contents of showFeaturesEdit as a double


% --- Executes during object creation, after setting all properties.
function showFeaturesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showFeaturesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadSecondImagePush.
function loadSecondImagePush_Callback(hObject, eventdata, handles)
% hObject    handle to loadSecondImagePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Function to get image
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.bmp;*.jpeg;*.png;*.gif','All Image Files';'*.*','All Files'}, 'Select Secong Image to Match');
fileName = fullfile(pathname, filename);

% read image
img2 = imread(fileName);

axes(handles.axes2)
imshow(rgb2gray(img2),[]);

handles.img2 = img2;
guidata(hObject, handles);



function matchThresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to matchThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of matchThresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of matchThresholdEdit as a double


% --- Executes during object creation, after setting all properties.
function matchThresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matchThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxRatioMatchEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxRatioMatchEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRatioMatchEdit as text
%        str2double(get(hObject,'String')) returns contents of maxRatioMatchEdit as a double


% --- Executes during object creation, after setting all properties.
function maxRatioMatchEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxRatioMatchEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadImagesToCalibratePush.
function loadImagesToCalibratePush_Callback(hObject, eventdata, handles)
% hObject    handle to loadImagesToCalibratePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%to get the dialog to select the directory of calibration images
path = uigetdir('cwd');

% to get the files in the directory selected
matFiles = dir2([path,'\*.jpg']);

%to get number of files in the directory
numFiles = length(matFiles);

%to sort file names accordingly
fileNames = natsortfiles({matFiles.name});

for nFile = 1 : numFiles
    fileName = strcat(path, '\', fileNames{nFile}); % to get filename
    
    imgNames{nFile} = fileName;

end
handles.imgNames = imgNames;
guidata(hObject, handles);
set(handles.processPush,'Enable','on'); % set process button vissibility on


% --- Executes on button press in loadImageFundamentalPush.
function loadImageFundamentalPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadImageFundamentalPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Function to get image
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.bmp;*.jpeg;*.png;*.gif','All Image Files';'*.*','All Files'}, 'Select Secong Image to Match');
fileName = fullfile(pathname, filename);

% read image
img2 = imread(fileName);

axes(handles.axes2)
imshow((img2),[]);

handles.img2 = img2;
guidata(hObject, handles);



function f1Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f1Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1Edit as text
%        str2double(get(hObject,'String')) returns contents of f1Edit as a double


% --- Executes during object creation, after setting all properties.
function f1Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f2Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f2Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f2Edit as text
%        str2double(get(hObject,'String')) returns contents of f2Edit as a double


% --- Executes during object creation, after setting all properties.
function f2Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f2Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f3Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f3Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f3Edit as text
%        str2double(get(hObject,'String')) returns contents of f3Edit as a double


% --- Executes during object creation, after setting all properties.
function f3Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f3Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f4Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f4Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f4Edit as text
%        str2double(get(hObject,'String')) returns contents of f4Edit as a double


% --- Executes during object creation, after setting all properties.
function f4Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f4Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f5Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f5Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f5Edit as text
%        str2double(get(hObject,'String')) returns contents of f5Edit as a double


% --- Executes during object creation, after setting all properties.
function f5Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f5Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f6Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f6Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f6Edit as text
%        str2double(get(hObject,'String')) returns contents of f6Edit as a double


% --- Executes during object creation, after setting all properties.
function f6Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f6Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f7Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f7Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f7Edit as text
%        str2double(get(hObject,'String')) returns contents of f7Edit as a double


% --- Executes during object creation, after setting all properties.
function f7Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f7Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f8Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f8Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f8Edit as text
%        str2double(get(hObject,'String')) returns contents of f8Edit as a double


% --- Executes during object creation, after setting all properties.
function f8Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f8Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f9Edit_Callback(hObject, eventdata, handles)
% hObject    handle to f9Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f9Edit as text
%        str2double(get(hObject,'String')) returns contents of f9Edit as a double


% --- Executes during object creation, after setting all properties.
function f9Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f9Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stitchPush.
function stitchPush_Callback(hObject, eventdata, handles)
% hObject    handle to stitchPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        path = uigetdir('cwd');
        files = dir2([path]);
        %to get number of files in the directory
        numFiles = length(files);
        
        %to sort file names accordingly
        fileNames = natsortfiles({files.name});
        
        for nFile = 1 : numFiles
            fileName = strcat(path, '\', fileNames{nFile}); % to get filename
            imgS{nFile} = imread(fileName);
            
        end
        set(handles.processPush, 'enable','on');
        handles.imgS = imgS;
        guidata(hObject, handles);


% --- Executes on button press in savePush.
function savePush_Callback(hObject, eventdata, handles)
% hObject    handle to savePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imgout = handles.img_out;
[a b] =uiputfile();
filename = strcat(b,a);
imwrite(imgout,filename);
set(handles.saveText,'visible','on');


% --- Executes on button press in resultCheck.
function resultCheck_Callback(hObject, eventdata, handles)
% hObject    handle to resultCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resultCheck


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over resultCheck.
function resultCheck_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to resultCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img_out = handles.img_out;
handles.img_out = img_out;
guidata(hObject,handles);


% --- Executes on button press in loadVideoPush.
function loadVideoPush_Callback(hObject, eventdata, handles)
% hObject    handle to loadVideoPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Function to get video
[filename, pathname] = uigetfile({'*.avi;*.mpeg','All Video Files';'*.*','All Files'}, 'Select a Video');
fileName = fullfile(pathname, filename);

h = waitbar(0,'Loading Video...');
       
    %create a video object
    vidObject = VideoReader(fileName);
    
    %get video frame height and width
    vidHeight = vidObject.Height;
    vidWidth = vidObject.Width;
    
    %take 200 frames in each video
    vidFrames = vidObject.NumberOfFrames;
    
    % Create a video as a 3 dimensional data
    video3Dmat = zeros(vidHeight, vidWidth, vidFrames);
    
    % Read one frame at a time using readFrame until 200 frames reached
    % Append data from each video frame to the mat array
    for nFrame = 1: vidFrames
        vid_image = read(vidObject, nFrame);
        img(:, :,:, nFrame) = vid_image ;
        waitbar(nFrame/vidFrames);
    end
    
    close(h)
    handles.img = img;
    guidata(hObject, handles);
    
    set(handles.processPush,'enable','on');
