function [feature_x, feature_y] = HarrisDetection(img, w, sigma, k, threshold)
    % This function implements Harris corner detection
	% Turn RGB image into gray scale image
	[row, col, dim] = size(img);
	I = rgb2gray(img);    
	I_double = double(I);
	R = zeros(row, col);

	% 1. Compute x and y derivatives of image
	IG = GaussianFunction(I_double, w, sigma);
	[I_x, I_y] = gradient(IG);

	% 2. Compute products of derivatives at every pixel
	I_x2 = I_x .^ 2;
	I_y2 = I_y .^ 2;
	I_xy = I_x .* I_y;

	% 3. Compute the sums of the products of derivatives at each pixel
	S_x2 = GaussianFunction(I_x2, w, sigma);
	S_y2 = GaussianFunction(I_y2, w, sigma);
	S_xy = GaussianFunction(I_xy, w, sigma);

	% 4. Define the matrix at each pixel
	for x = 1 : row
		for y = 1 : col
			M = [S_x2(x, y) S_xy(x, y); S_xy(x, y) S_y2(x, y)];

			% 5. Compute the reponse of the detector at each pixel
			R(x, y) = det(M) - k * (trace(M) ^ 2);
		end
    end
    
	% 6. Threshold on value of R; compute nonmax supression
	R_thres = (R > threshold);

	mask = [1 1 1; 1 0 1; 1 1 1];
	R_NMS = (R > imdilate(R, mask));
	R_LM = R_thres & R_NMS;
	[feature_x_result, feature_y_result] = find(R_LM);
	DrawPoint(img, feature_x_result, feature_y_result);
	% imshow(R_LM);

	% Remove low contrast
    D = GaussianFunction(I_double, w, 2 * sigma) - GaussianFunction(I_double, w, sigma);
    D = D + 0.5 * (gradient(D));
    disp(D);
    % Remove edge
	% Parameters of derivative function
    % f(1, 0) - 2f(0, 0) + f(-1, 0)
	x2 = [1 -2 1];
    % f(0, 1) - 2f(0, 0) + f(0, -1)
	y2 = [1; -2; 1];
    % f(-1, -1) - f(-1, 1) - f(1, -1) + f(1, 1) / 4
	xy = [1 0 -1; 0 0 0; -1 0 1];

	% Thresholds
	edge_threshold = ((10 + 1) ^ 2) / 10;
	contrast = abs(filter2(fspecial('gaussian', 5), I_double) - I_double);
	contrast_threshold = 0.03;

	feature_x = [];
	feature_y = [];

	for i = 1:numel(feature_x_result)
		x = feature_x_result(i);
		y = feature_y_result(i);

		% Remove boundary
		if ((x > 7) && (x <= (col - 7)) && (y > 7) && (y <= col - 7))
			% Remove edge
			D_x2 = sum(I_double(y, x - 1:x + 1) .* x2);
			D_y2 = sum(I_double(y - 1:y + 1, x) .* y2);
			D_xy = sum(sum(I_double(y - 1:y + 1, x - 1:x + 1) .* xy)) / 4;

			Tr_Hessian = D_x2 + D_y2;
			Det_Hessian = D_x2 * D_y2 - D_xy ^ 2;

			ratio = (Tr_Hessian ^ 2) / Det_Hessian;

			if ((Det_Hessian >= 0) && (ratio < edge_threshold) && (D(y, x) < contrast_threshold))
				feature_x = [feature_x; x];
				feature_y = [feature_y; y];
            end
        end
    end
    DrawPoint(img, feature_x, feature_y);
end