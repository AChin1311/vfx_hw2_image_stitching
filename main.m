function main()
    read_cache = 0;
    save_cache = 0;
    
    % Load Images
	disp('Loading Images');
    image_serial = 'grail/';
    directory = ['image/' image_serial];
	output_filename = [image_serial '_stitched.png'];
    files = dir([directory, '*.jpg']);
    img_num = length(files);
    img_info = imfinfo([directory, files(1).name]);
    imgy = img_info.Height;
    imgx = img_info.Width;
    img_array = {};
    
    % Load focal length
    fileID = fopen([directory 'focal_len.txt'], 'r');
    focals = fscanf(fileID, '%f');
    fclose(fileID);
    
    % Features detection
    disp('Features detection');

    for i = 1 : img_num
        ImagePath = [directory, files(i).name];
        img = imread(ImagePath); 
        warpimg{i} = warpFunction(img, focals(i));
        %img_array(:, :, :, i) = warpimg{i};
    end

    % Features detection
    disp('Features detection');
    for i = 1:5
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
            [pos, orient, desc] = SIFTdescriptor(img, fx, fy);
            disp(size(pos));

            
            if (save_cache)
                save(sprintf('image/%s/mat/fx_%02d.mat', image_serial, i), 'fx');
                save(sprintf('image/%s/mat/fy_%02d.mat', image_serial, i), 'fy');
                save(sprintf('image/%s/mat/pos_%02d.mat', image_serial, i), 'pos');
                save(sprintf('image/%s/mat/orient_%02d.mat', image_serial, i), 'orient');
                save(sprintf('image/%s/mat/desc_%02d.mat', image_serial, i), 'desc');
            end
        end
        % DrawArrow(img, pos(:, 1), pos(:, 2), orient);
        
        fxs{i} = fx;
        fys{i} = fy;
        poss{i} = pos;
        orients{i} = orient;
        descs{i} = desc;
    end
    
    % DrawPoint(warpimg{1}, fys{1}, fxs{1});
    % DrawPoint(warpimg{2}, fys{2}, fxs{2});

    for i = 1 : 4
        match = ransac(warpimg{i}, descs{i}, poss{i}, warpimg{i + 1}, descs{i + 1}, poss{i + 1});
        trans = matchImage(match, poss{i}, poss{i + 1});
        blendImage(warpimg{i}, warpimg{i + 1}, trans);
    end
    
    % DrawPoint(warpimg{1}, poss{1}(match(:, 1), 2), poss{1}(match(:, 1), 1));
    % DrawPoint(warpimg{2}, poss{2}(match(:, 2), 2), poss{2}(match(:, 2), 1));
    % Features matching
    disp('Features matching');

    
    % Images matching
    disp('Images matching');

    % Images blending
    disp('Images blending');
end