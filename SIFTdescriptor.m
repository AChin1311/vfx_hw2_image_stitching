function [pos, orient, desc] = SIFTdescriptor(img, feature_x, feature_y)

    orient = [];
    pos = [];
    desc = [];
    
    % Orientation assignment
    I = double(rgb2gray(img));
    L = filter2(fspecial('gaussian', [5 5]), I);
    [y_len x_len dim] = size(img);
    
    % Gradient magnitude
    % Dx = L(x+1, y) - L(x-1, y); Dy = L(x, y+1) - L(x, y-1)
    % m(x, y) = sqrt(Dx^2 - Dy^2)
    Dx = 0.5*(L(2:(end-1), 3:end)-L(2:(end-1), 1:(end-2)));
    Dy = 0.5*(L(3:end, 2:(end-1))-L(1:(end-2), 2:(end-1)));
    magnitude = zeros(y_len, x_len);
    magnitude(2:(end-1), 2:(end-1)) = sqrt(Dx.^2 + Dy.^2);
    
    % Gradient orientation
    gradient = zeros(y_len, x_len);
    gradient(2:(end-1), 2:(end-1)) = atan2(Dy, Dx);
    gradient(find(gradient == pi)) = -pi;
    
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
            add = sum(weighted_mag.*max(1-abs(diff)/step,0));
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
        % find the maxpeak
        [maxPeakValue peakIndex] = max(peaks);
        peakValue = maxPeakValue;
        while(peakValue > maxPeakValue*0.8)
            max_orient = orient_ref(peakIndex)+pi/36;
            pos = [pos; [x, y]];
            orient = [orient; max_orient];
            % find the next peak
            peaks(peakIndex) = 0;
            [peakValue peakIndex] = max(peaks);
        end
    end

    disp(size(pos));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spin_step = pi/4;
    orient_angles = [-pi:spin_step:(pi-spin_step)];
    
    grid_step = 4;
    feat_window = 2*grid_step;
    
    [cor_x cor_y] = meshgrid( [-6:grid_step:6] );
    feat_grid = [cor_x(:) cor_y(:)]';
    
    [cor_x cor_y] = meshgrid( [-(2*grid_step-0.5):(2*grid_step-0.5)] );
    feat_samples = [cor_x(:) cor_y(:)]';
    
    for k = 1:size(orient)
        x = pos(k,1);
        y = pos(k,2); 
        
        % rotate the grid
        M = [cos(orient(k)) -sin(orient(k)); sin(orient(k)) cos(orient(k))];
        rot_grid = M*feat_grid + repmat([x; y], 1, 16);
        rot_samples = M*feat_samples + repmat([x; y], 1, 256);
        feat_desc = zeros(1,128);
        
        % sampled over 16x16 array of locations
        start = 0;
        for s = 1:4
            for t = 1:4
                block_bins = zeros(1, 8);
                for i = 1:4
                    start = (s-1)*64+(t-1)*4+1+(i-1)*16;
                    pixel_x = round(rot_samples(1, start:(start+3)));
                    pixel_y = round(rot_samples(2, start:(start+3)));
                    b_x = zeros(1, 6);
                    b_y = zeros(1, 6);
                    b_x(1) = 2*pixel_x(1) - pixel_x(2);
                    b_x(2:5) = pixel_x(1:4);
                    b_x(6) = 2*pixel_x(4) - pixel_x(3);
                    b_y(1) = 2*pixel_y(1) - pixel_y(2);
                    b_y(2:5) = pixel_y(1:4);
                    b_y(6) = 2*pixel_y(4) - pixel_y(3);
%                     disp(pixel_x);
%                     disp(pixel_y);
%                     disp(b_x);
%                     disp(b_y);
                    
                    if isempty(find(b_x <= 1)) && isempty(find(b_y <= 1)) && ...
                       isempty(find(b_x >= x_len-1)) && isempty(find(b_y >= y_len-1))
                        %disp('ok');
                    else     
                        %disp('skip');
                        if i == 4
                            disp('flush');
                            samples_desc((s-1)*4+t, 1:8) = [1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8];
                        end
                        continue;
                    end
                    
                    bin_ref = [-pi:pi/4:(pi-pi/4)];
                    for d = 2:5
                        px_Dx = 0.5*(L(b_y(d), b_x(d+1)) - L(b_y(d), b_x(d-1)));
                        px_Dy = 0.5*(L(b_y(d+1), b_x(d)) - L(b_y(d-1), b_x(d)));
                        px_mag = sqrt(px_Dx^2 + px_Dy^2);
                        px_grad = atan(px_Dy/px_Dx);
                        disp([px_Dy/px_Dx, px_Dx^2 + px_Dy^2]);
                        if px_grad == pi
                            px_grad = -pi;
                        end
                        if isnan(px_grad)
                            continue;
                        end
                        %disp(px_mag);
                        %disp(px_grad);

                        b = 1;
                        while px_grad > bin_ref(b)
                            b = b+1;
                        end
                        %disp(b);
                        block_bins(b-1) = block_bins(b-1)+px_mag; 
                    end
                    if(i == 4)
                        %disp('flush');
                        samples_desc((s-1)*4+t, 1:8) = block_bins/norm(block_bins);
                    end
                end
            end
        end
        disp(samples_desc);
        % clip values larger than 0.2
        samples_desc(find(samples_desc > 0.2)) = 0.2;
        % renormalize
        samples_desc = samples_desc/norm(samples_desc);
        desc = [desc; reshape(samples_desc, 1, 128)];

    end
end