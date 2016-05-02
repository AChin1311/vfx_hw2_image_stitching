function [orient, pos, desc] = SIFTdescriptor(img, feature_x, feature_y)

    orient = [];
    pos = [];
    desc = [];

    I = rgb2gray(img);
    % Orientation assignment
    L = filter2(fspecial('gaussian', [5 5]), I);
    [y_len x_len] = size(I);
    
    % Gradient magnitude
    % Dx = L(x+1, y) - L(x-1, y); Dy = L(x, y+1) - L(x, y-1)
    % m(x, y) = sqrt(Dx^2 - Dy^2)
    Dx = 0.5*L(2:(end-1), 3:end)-L(2:(end-1), 1:(end-2));
    Dy = 0.5*L(3:end, 2:(end-1))-L(1:(end-2), 2:(end-1));
    magnitude = zeros(y_len, x_len);
    magnitude(2:(end-1), 2:(end-1)) = sqrt(Dx.^2 + Dy.^2);
    
    % Gradient orientation
    grad = zeros(y_len, x_len);
    grad(2:(end-1), 2:(end-1)) = atan2(Dy, Dx);
    grad(find(grad == pi)) = -pi;
    
    % set up histogram
    step = pi/18;
    orient_ref = zeros(36, 1);
    for i = 1:36
        orient_ref(i) = -pi+(i-1)*step; 
    end
    
    % Weighted by 2D gaussian kernel
    % sigma = 1.5 * scale of the keypoint
    winSize = 7;
    half = floor(winSize/2);
    mask = fspecial('gaussian', [winSize winSize], 1.5*0.5);
    
    % Time to vote!
    for k = 1:numel(feature_x)
        x = feature_x(k);
        y = feature_y(k);
        
        weighted_mag = mask .* magnitude((y-half):(y+half), (x-half):(x+half));
        window_orient = L((y-half):(y+half), (x-half):(x+half));
        histogram = zeros(36, 1);
        for bin = 1:36
            diff = mod(window_orient-orient_ref(bin)+pi, 2*pi) - pi;
            add = sum(weighted_mag.*max(1-abs(diff)/step,0))
            histogram(bin) = histogram(bin)+sum(add);
        end
        
        peaks = zeros(36, 1);
        for i = 2:35
            if(histogram(i) > histogram(i-1) && histogram(i) > histogram(i+1))
                peaks(i) = histogram(i);
            end
        end
        if(histogram(1) > histogram(36) && histogram(1) > histogram(2))
            peaks(1) = histogram(1);
        end
        if(histogram(36) > histogram(35) && histogram(36) > histogram(1))
            peaks(36) = histogram(36);
        end
        
        % find maxpeak
        [maxPeakValue peakIndex] = max(peaks);
%         disp('maxPeakValues');
%         disp(maxPeakValue);
%         disp('peakIndex');
%         disp(peakIndex);
        peakValue = maxPeakValue;
        
        while(peakValue > maxPeakValue*0.8)
            A = [];
            b = [];
            for j = -1:1
                A = [A; (orient_ref(peakIndex)+step*j).^2 (orient_ref(peakIndex)+step*j) 1];
                bin = mod(peakIndex + j + 36 - 1, 36) + 1;
                b = [b; peaks(bin)];
            end
            c = pinv(A)*b;
            max_orient = -c(2)/(2*c(1));
            while(max_orient < -pi)
                max_orient = max_orient + 2*pi;
            end
            while(max_orient >= pi)
                max_orient = max_orient - 2*pi;
            end

            pos = [pos; [x, y]];
            orient = [orient; max_orient];

            % find the next peak
            peaks(peakIndex) = 0;
            [peakValue peakIndex] = max(peaks);
        end
    end
%     disp('pos');
%     disp(pos);
%     disp('ori');
%     disp(orient);
%     disp(size(pos));
%     disp(size(orient));
    
    
    theta = pi/4;
    orient_angles = [-pi:theta:(pi-theta)];
    
    [x_coords y_coords] = meshgrid( [-6:4:6] );
    grid = [x_coords(:) y_coords(:)]';
    [x_coords y_coords] = meshgrid( [-(2*4-0.5):(2*4-0.5)] );
    feat_samples = [x_coords(:) y_coords(:)]';

%     cor = [-6, -2, 2, 6];
%     for i = 1:4
%         for j = 1:4
%             grid(1, (i-1)*4+j) = cor(i);
%             grid(2, (i-1)*4+j) = cor(j);
%         end
%     end
%     
%     sample_cor = [-7.5:7.5];
%     for i = 1:15
%         for j = 1:15
%             feat_samples(1, (i-1)*15+j) = sample_cor(i);
%             feat_samples(2, (i-1)*15+j) = sample_cor(j);
%         end
%     end
    feat_window = 8;  
    
    for k = 1:size(orient)
        x = pos(k,1);
        y = pos(k,2);
        % Rotate grid
        M = [cos(orient(k)) -sin(orient(k)); sin(orient(k)) cos(orient(k))];
        rot_grid = M*grid + repmat([x; y],1,size(grid,2));
        rot_samples = M*feat_samples + repmat([x; y],1,size(feat_samples,2));
        
%         M = [cos(orient(k)) -sin(orient(k)); sin(orient(k)) cos(orient(k))];
%         for i = 1:16
%                 tmp1(1, i) = x;
%                 tmp1(2, i) = y;
%         end
%         rot_grid = M*grid + tmp1;
%         for i = 1:225
%                 tmp2(1, i) = x;
%                 tmp2(2, i) = y;
%         end
%         rot_samples = M*feat_samples + tmp2;
   
        feat_desc = zeros(1,128);
        for s = 1:size(rot_samples,2)
            x_sample = rot_samples(1,s);
            y_sample = rot_samples(2,s);
            for i = 1:3
                X(i, 1) = x_sample-1;
                X(i, 2) = x_sample;
                X(i, 3) = x_sample+1;
                Y(1, i) = y_sample-1; 
                Y(2, i) = y_sample;
                Y(3, i) = y_sample+1;
            end
            
            G = interp2(I, X, Y, '*linear');
            G(find(isnan(G))) = 0;
            
            Dx = 0.5*(G(2,3) - G(2,1));
            Dy = 0.5*(G(3,2) - G(1,2));
            mag_sample = sqrt(Dx^2 + Dy^2);
            grad_sample = atan2(Dy, Dx);
            if grad_sample == pi
                grad_sample = -pi;
            end
            
            x_weight = max(1 - (abs(rot_grid(1,:) - x_sample)/4), 0);
            y_weight = max(1 - (abs(rot_grid(2,:) - y_sample)/4), 0); 
            for i = 1:8
                tmp3(i, :) = x_weight.*y_weight;
            end
            pos_weight = reshape(tmp3,1,128);
            
            diff = mod(grad_sample - orient(k) - orient_angles + pi, 2*pi) - pi;
            orient_weight = max(1 - abs(diff)/theta,0);
            tmp4 = [];
            for i =1:16
                tmp4 = [tmp4 orient_weight];
            end
            orient_weight = tmp4;         
            g = exp(-((x_sample-x)^2+(y_sample-y)^2)/(2*feat_window^2))/(2*pi*feat_window^2);
            feat_desc = feat_desc + pos_weight.*orient_weight*g*mag_sample;
        end
        feat_desc = feat_desc / norm(feat_desc);

        feat_desc(find(feat_desc > 0.2)) = 0.2;
        feat_desc = feat_desc / norm(feat_desc);

        desc = [desc; feat_desc];
    end
    disp(desc);
end