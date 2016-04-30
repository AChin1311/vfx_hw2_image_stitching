function imout = blendImage(img1, img2, trans)
	[row1, col1, dim1] = size(img1);
	[row2, col2, dim2] = size(img2);
	imout = zeros(row1 + abs(trans(2)), col1 + abs(trans(1)), dim1);
end