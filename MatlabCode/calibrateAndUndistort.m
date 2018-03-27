function calibrateAndUndistort(  )
%this function is used to calibrate camera using calibration pattern 
%images and undistort the image captured by camera

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

%detect checkboard ppints
[imagePoints,boardSize,~] = detectCheckerboardPoints(imgNames);

% Generate world coordinates of the corners of the squares. Square size is in millimeters
squareSize = 29;
worldPoints = generateCheckerboardPoints(boardSize,squareSize);

% Calibrate the camera
cameraParams = estimateCameraParameters(imagePoints,worldPoints);

figure; showExtrinsics(cameraParams);
                              
% Remove lens distortion and display results.
I = imread(imgNames{10});
J1 = undistortImage(I,cameraParams);

figure; imshowpair(I,J1,'montage');
title('Original Image (left) vs. Corrected Image (right)');
                              
                              