function main()
    read_cache = 0;
    save_cache = 1;

    % Load Images
	disp('Loading Images');
    image_serial = 'book2';
    directory = ['image/' image_serial '/'];
	output_filename = [image_serial '_stitched.png'];
    files = dir([directory, '/*.jpg']);
    img_num = length(files);
    
    % Load focal length
    disp('Loading focal length file')
    fileID = fopen([directory 'focal_len.txt'], 'r');
    focals = fscanf(fileID, '%f');
    fclose(fileID);
    
    % Warping images
    disp('Warping Images')
    for i = 1 :img_num
        ImagePath = [directory, files(i).name];
        img = imread(ImagePath); 
        disp(i);
        warpimg{i} = warpFunction(img, focals(i));
    end

    % Features detection
    disp('Features detection');
    for i = 1:img_num
        if read_cache
                load(sprintf('image/%s/mat/fx_%02d.mat', image_serial, i));
                load(sprintf('image/%s/mat/fy_%02d.mat', image_serial, i));
                load(sprintf('image/%s/mat/orient_%02d.mat', image_serial, i));
                load(sprintf('image/%s/mat/pos_%02d.mat', image_serial, i));
                load(sprintf('image/%s/mat/desc_%02d.mat', image_serial, i));
        else
            img = warpimg{i};
            [fx, fy] = HarrisDetection(img, 5, 1, 0.04, 3);   
            disp('key point');
            disp(size(fx));
            %[pos, orient, desc] = SIFTdescriptor(img, fx, fy);
            [pos, orient, desc] = WTFdescriptor(img, fx, fy);
            disp(size(pos));
            
            if (save_cache)
                save(sprintf('image/%s/mat/fx_%02d.mat', image_serial, i), 'fx');
                save(sprintf('image/%s/mat/fy_%02d.mat', image_serial, i), 'fy');
                save(sprintf('image/%s/mat/pos_%02d.mat', image_serial, i), 'pos');
                save(sprintf('image/%s/mat/orient_%02d.mat', image_serial, i), 'orient');
                save(sprintf('image/%s/mat/desc_%02d.mat', image_serial, i), 'desc');
            end
        end
        
        fxs{i} = fx;
        fys{i} = fy;
        poss{i} = pos;
        orients{i} = orient;
        descs{i} = desc;
    end
    
    % Features matching
    disp('Features matching');
    for i = 1 : img_num - 1
        match{i} = ransac(descs{i}, poss{i}, descs{i + 1}, poss{i + 1});
    end
    
    % Images matching
    disp('Images matching');
    for i = 1 : img_num - 1
        trans{i} = matchImage(match{i}, poss{i}, poss{i + 1});
    end
    
    % Images blending
    disp('Images blending');
    imout = warpimg{1};
    for i = 2 : img_num
        imout = blendImage(imout, warpimg{i}, trans{i - 1});
    end

    % Saving panorama
    disp('Saving panorama');
    imshow(imout);
    imwrite(imout, output_filename);
end