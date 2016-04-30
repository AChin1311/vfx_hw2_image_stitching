function main()
    % Load Images
	disp('Loading Images');
    image_serial = 'parrington/';
    directory = ['image/' image_serial];
	output_filename = [image_serial '_stitched.png'];
    files = dir([directory, '*.jpg']);
    img_num = length(files);
    img_info = imfinfo([directory, files(1).name]);
    imgy = img_info.Height;
    imgx = img_info.Width;
    img_array = zeros(imgy, imgx, 3, img_num);
    
    % Load focal length
    fileID = fopen([directory 'focal_len.txt'], 'r');
    focals = fscanf(fileID, '%f');
    %disp(focals);
    fclose(fileID);
    
    for i = 1 : img_num
        ImagePath = [directory, files(i).name];
        img = imread(ImagePath); 
        warpimg = warpFunction(img, focals(i));
        img_array(:, :, :, i) = warpimg;
    end
    
    % Features detection
    disp('Features detection');
    %testing
    img = warpimg %last one
    [feature_x, feature_y] = HarrisDetection(img, 5, 1, 0.04, 3);
    disp(feature_x);
    SIFTdescriptor(img, feature_x, feature_y);
    
    % Features matching
    disp('Features matching');

    % Images matching
    disp('Images matching');

    % Images blending
    disp('Images blending');
    
end