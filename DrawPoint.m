function DrawPoint(img, x, y)
    disp('Draw Point');
    figure('Position', [50 50 size(img,2) size(img,1)])
    image(img);
    hold on;
    plot(y(:), x(:), 'r.', 'MarkerSize', 20);
    hold off;
end