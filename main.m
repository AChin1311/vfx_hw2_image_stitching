function main()

	image_serial = 'prtn';
	directory = ['images/' image_serial];
	output_filename = [image_serial '_stitched.png'];
    
    % Images loading
    disp('Images loading');

    % Features detection
    disp('Features detection');

    % test of warp function
    img = imread('image/parrington/prtn00.jpg');
    warpFunction(img, 704.916);

    % Features matching
    disp('Features matching');

    % Images matching
    disp('Images matching');

    % Images blending
    disp('Images blending');
    
end