function main()
    % Load Images
	disp('Loading Images');
    image_serial = 'parrington/';
    directory = ['image/' image_serial];
	output_filename = [image_serial '_stitched.png'];
    files = dir([directory, '*.jpg']);
    img_num = length(files);
    img_info = imfinfo([directory, files(1).name]);
    imgRow = img_info.Height;
    imgCol = img_info.Width;
    disp([img_num imgRow imgCol]);
    img_array = zeros(imgRow, imgCol, 3, img_num);
    
    % Load focal length
    fileID = fopen([directory 'focal_len.txt'], 'r');
    focals = fscanf(fileID, '%f');
    disp(focals);
    fclose(fileID);
    
    for i = 1 : img_num
        ImagePath = [directory, files(i).name];
        img = imread(ImagePath);    
        warpimg = warpFunction(img, focals(i));
        disp(focals(i));
        img_array(:, :, :, i) = warpimg;
    end
    
    % Features detection
    disp('Features detection');

    
    % Features matching
    disp('Features matching');

    % Images matching
    disp('Images matching');

    % Images blending
    disp('Images blending');
    
end