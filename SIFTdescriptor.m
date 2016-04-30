function SIFTdescriptor(img, feature_x, feature_y)
    % Orientation assignment
    disp('SIFTdescriptor');
    sigma = 1;
    L = GaussianFunction(rgb2gray(img), 5, sigma);
    [imgy imgx dim] = size(img);
    
    % Gradient magnitude
    % Dx = L(x+1, y) - L(x-1, y); Dy = L(x, y+1) - L(x, y-1)
    % m(x, y) = sqrt(Dx^2 - Dy^2)
    Dx = L(2:(end-1), 3:end) - L(2:(end-1), 1:(end-2));
    Dy = L(3:end, 2:(end-1)) - L(1:(end-2), 2:(end-1));
    magnitude = zeros(imgy, imgx);
    magnitude(2:(end-1), 2:(end-1)) = sqrt(Dx.^2 + Dy.^2);
    
    % Gradient orientation
    orient = zeros(imgy, imgx);
    orient(2:(end-1), 2:(end-1)) = round(atan2(Dy, Dx)*(180/pi)+180);
    
    
    % Weighted by 2D gaussian kernel
    % sigma = 1.5 * scale of the keypoint
    winSize = 7;
    half = floor(winSize/2);
    mask = fspecial('gaussian', [winSize winSize], 1.5*sigma);
    
    % Time to vote!
    for k = 1:size(feature_x)
        x = feature_x(k);
        y = feature_y(k);
        disp([y x]);
        
        weighted_magnitude = mask .* magnitude((y-half):(y+half), (x-half):(x+half));
        window_orient = orient((y-half):(y+half), (x-half):(x+half));
        buckets = zeros(36, 1);
        for w_y = 1:winSize
            for w_x = 1:winSize
                b = min(floor(window_orient(w_y, w_x)/10)+1, 36);
                
                buckets(b) = buckets(b) + window_orient(w_y, w_x)*weighted_magnitude(w_y, w_x);
            end
        end
        disp(buckets);
        % Check the voting result
        for i = 2:35
            if buckets(i) < buckets(i-1) || buckets(i) < buckets(i+1)
                buckets(i) = 0;
            end
        end
        % handle the first and last bin
        if buckets(1) < buckets(36)
            buckets(1) = 0;
        end
        if buckets(36) < buckets(1)
            buckets(36) = 0;
        end
        % find the peaks
        [maxPeakValue peakIndex] = max(buckets);
        peakValue = maxPeakValue;
        while(peakValue > maxPeakValue*0.8)
            orientation = -pi + (peakIndex-1)*pi/18;
            A = [(orientation-pi/18)^2 (orientation-pi/18) 1;...
                  orientation.^2        orientation        1;...
                 (orientation+pi/18)^2 (orientation+pi/18) 1];
            bin = [mod(peakIndex+34, 36) + 1; mod(peakIndex+35, 36) + 1; mod(peakIndex+36, 36) + 1];
            b = [buckets(bin(1)) buckets(bin(2)) buckets(bin(3))];
            c = pinv(A)*b;
            des_orient = -c(2)/(2*c(1));
            if (des_orient < -pi)
                n = floor((des_orient+pi)/(2*pi))+1;
                des_orient = des_orient + 2*n*pi;
            end
            if (des_orient >= pi)
                n = floor((des_orient-pi)/(2*pi))+1;
                des_orient = des_orient - 2*n*pi;
            end
            orient = [orient; des_orient];
            pos = [pos; [x, y]];
            
            % find the next peak
            [peakValue peakIndex] = max(Buckets);
        end
    end
    % Local image descriptor
    
end