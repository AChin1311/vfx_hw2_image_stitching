function imout = GaussianFunction(img, w, sigma)
	% Gaussian function
	f = fspecial('gaussian', [w w], sigma);
	imout = filter2(f, img);
end