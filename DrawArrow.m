function DrawArrow(img, x, y, ori)
    angle = round(ori(:) ./pi .*180);
    angle(find(angle == -180)) = 180;
    disp(angle);
    for i = 1:size(angle);
        text(y(i), x(i), '\leftarrow','rotation', angle(i));
    end
end
