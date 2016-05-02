function imout = blendImage(im1, im2, trans)
	%[row1, col1, dim1] = size(img1);
	%[row2, col2, dim2] = size(img2);
	%imout = zeros(row1 + abs(trans(2)), col1 + abs(trans(1)), dim1);

	[row1, col1, channel] = size(im1);
    [row2, col2, channel] = size(im2);
    imout = zeros(row1+trans(2), col1+abs(trans(1)), channel);

    blendWidth = trans(1);

    % merge by 'plus'
    for y = 1:row1
        for x = 1:col1
            imout(y,x,:) = im1(y,x,:);
        end
    end
    %for y = 1:row2
     %   for x = 1:col2 
      %      x1 = x + size(imout, 2) - col2;
       %     imout(y,x1,:) = im2(y,x,:);
        %end
    %end

    imshow(imout, 'DisplayRange', [0, 255]);
end