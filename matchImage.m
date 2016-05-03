function trans = matchImage(match, pos1, pos2)
	% Solve 2D translation
	n = size(match, 1);
	A = zeros(n * 2 + 1, 3);
	b = zeros(n * 2 + 1, 1);

	for i = 1 : n
		% Fill x and x' into matrices
		A(i, 1) = 1;
		A(i, 3) = pos1(match(i, 1), 1);
		b(i, 1) = pos2(match(i, 2), 1);

		% Fill y and y' into matrices
		A(i + n, 2) = 1;
		A(i + n, 3) = pos1(match(i, 1), 2);
		b(i + n, 1) = pos2(match(i, 2), 2);
	end
	A(n * 2 + 1, 3) = 1;
	b(n * 2 + 1, 1) = 1;
    
    disp('match');
    disp(match);
    disp('pos');
    disp(pos1);

	% Solve transaltion
	result = round(A \ b);
	trans = result(1 : 2);
    disp('result');
    disp(result);
    disp(trans);
end