function [feature_x, feature_y, R] = HarrisDetection(img, w, sigma, k, threshold)
    % This function implements Harris corner detection

	% Turn RGB image into gray scale image
	[row, col, dim] = size(img);
    disp([row col]);
	I = rgb2gray(img);
	I = double(I);
	R = zeros(row, col);

	% 1. Compute x and y derivatives of image
	IG = GaussianFunction(I, 5, 1);
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

	[feature_x, feature_y] = find(R_lm);
end