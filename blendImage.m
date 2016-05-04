function imout = blendImage(img1, img2, trans)
	[row1, col1, channel] = size(img1);
	[row2, col2, channel] = size(img2);
	imout = zeros(row1 + abs(trans(2)), col1 + abs(trans(1)), channel, 'uint8');
    
    blendWidth = col2 - trans(1);  
    
    % Alpha layer
    alpha1 = ones(row1, col1); 
    alpha2 = ones(row2, col2);
    
    for i = 1 : row2
        alpha2(i, (end - blendWidth):end) = 1:(-1 / blendWidth):0;
        for j = 1 : col2
            imout(i, j, :) = img2(i, j, :) * alpha2(i, j);
        end
    end
    
    for i = 1 : row1
        alpha1(i, 1:blendWidth + 1) = 0:(1 / blendWidth):1;
        for j = 1 : col1
            tmp = imout(i + abs(trans(2)), j + abs(trans(1)), :);
            imout(i + abs(trans(2)), j + abs(trans(1)), :) = tmp + img1(i, j, :) * alpha1(i, j);
        end
    end
    imshow(imout);
end