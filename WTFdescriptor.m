function [pos, orient, desc] = WTFdescriptor(img, fx, fy)
    pos = [];
    orient = [];
    desc = [];
    for s = 1:size(fx)
        x = fx(s);
        y = fy(s);
        pos = [pos; [x, y]];
        
        feat_desc = zeros(49, 3);
        k = 1;
        for i = -3:3
            for j = -3:3
                feat_desc(k, :) = [img(y+i, x+j, 1), img(y+i, x+j, 2), img(y+i, x+j, 3)];
                k = k +1;
                %git disp(feat_desc(k-1, :));
            end
        end
        feat_desc = reshape(feat_desc, 1, []);
        desc = [desc; feat_desc];
    end
end