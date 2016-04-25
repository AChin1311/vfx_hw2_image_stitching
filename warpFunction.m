function warpFunction(img, focal_length)
	% This function project the original image on a cylinder

	imout = zeros(size(img), 'uint8');
	mid = size(img) / 2;
	% Let s = focal length gives less distortion
	s = focal_length;

	for y = 1:size(img, 1)
		for x = 1:size(img, 2)
			x_dis = x - mid(2);
			y_dis = y - mid(1);

			theta = atan(x_dis / focal_length);
			h = y_dis / sqrt(x_dis ^ 2 + focal_length ^ 2);

			x_p = round(s * theta) + mid(2);
			y_p = round(s * h) + mid(1);

			imout(y_p, x_p, :) = img(y, x, :);
		end
	end

	imshow(imout);
end