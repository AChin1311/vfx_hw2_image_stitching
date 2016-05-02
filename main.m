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
    
    % Features detection
    disp('Features detection');
    for i = 1 : 2
        ImagePath = [directory, files(i).name];
        img = imread(ImagePath); 
        warpimg{i} = warpFunction(img, focals(i));
        img_array(:, :, :, i) = warpimg{i};
        [feature_x, feature_y] = HarrisDetection(warpimg{i}, 5, 1, 0.04, 3);
        [orient{i}, pos{i}, desc{i}] = SIFTdescriptor(warpimg{i}, feature_x, feature_y);
    end
    
    match = ransac(desc{1}, pos{1}, desc{2}, pos{2});
    trans = matchImage(match, pos{1}, pos{2});
    
    %testing
    %img = warpimg; %last one
    %[feature_x, feature_y] = HarrisDetection(img, 5, 1, 0.04, 3);
    %[orient, pos, desc] = SIFTdescriptor(img, feature_x, feature_y);
    %disp(pos);
    
    % Features matching
    disp('Features matching');

    % Images matching
    disp('Images matching');

    % Images blending
    disp('Images blending');
    blendImage(warpimg{1}, warpimg{2}, trans);
    
end