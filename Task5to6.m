clear; close all;
% Task 5: Robust method --------------------------
%Author: Jiawei Liu 
%18741916@students.lincoln.ac.uk

%----------------------IMPORTANT--------------------------------
% before run this program:
% 1.please create a file named"output1" to store all binary image;
% 2.plese create a file named"output2" to store result image;
% 3.please create a excel file named'data.xlsx' to store all evulation score.
% 4.Put this"Task5to6.m" file and the image that need to be processed under a file directory
% 5.all result will be print in command window, please check it .
% 
% best wish,
% Jiawei Liu
%--------------------  IMPORTANT--------------------------------------

data = [];
for k = 1 : 15 % use a For loop that get all images value(15 images)
    % get image file name by k
    filename = fullfile(pwd, sprintf('IMG_%02d.png', k));
    filename1 = fullfile(pwd, 'output1', sprintf('IMG_%02d.png', k));
    filename2 = fullfile(pwd, 'output2', sprintf('IMG_%02d.png', k));
    
    I = imread(filename);
    figure; imshow(I, []); title('origin image');
    % gray image
    I_grayscale = rgb2gray(I);
    figure; imshow(I_grayscale, []); title('gray image');
    % enhance
    I_rescale = imresize(I_grayscale,[512 NaN]);
    J = 255*im2double(I_rescale);
    mi = min(min(J));
    ma = max(max(J));
    I_Enhance = imadjust(I_rescale,[mi/255; ma/255],[0 1]);
    figure; imshow(I_Enhance, []); title('enhance image');
    % segmet 1: The cells are segment, all details are in the sub-function,
    % please pull down to check.
    res1 = main_seg1(I);
    % segmet 2: bacteria are segment,all details are in the sub-function
    %please pull down to check.
    res2 = main_seg2(I,res1);
    res1 = imresize(res1, [size(I,1) size(I,2)], 'bilinear');
    res2 = imresize(res2, [size(I,1) size(I,2)], 'bilinear');
    res = logical(res1 + res2);
    % show result image that have been sgemented
    figure; imshow(res, []); title('binary image');
    imwrite(res, filename1);
    % get result
    res = zeros(size(res1));
    % set label
    res(res1) = 1;
    res(res2) = 2;
    % result = label2rgb(res,@jet,[.0 .0 .0]);
    outputs = res;
    flag = 0;
    % If  only have background and bacteria, you need to change the flag to 0, 1
    if max(unique(outputs(:))) == 2 && length(unique(outputs(:)))==2
        flag = 1;
        outputs(outputs==2) = 1;
    end% give red blood cell and bacteria to a color
    L_res = label2rgb(res(:,:,1), 'prism','k','shuffle');
    imwrite(L_res, filename2);
    % dispaly result and save to image file
    figure; imshow(L_res, []); title('segment result');
    print(gcf, '-djpeg', '-r200', 'final image');

%---------------------------TASK 6-----------------------------------------

    % Task 6: Performance evaluation -----------------
    % Step 1: Load ground truth data
    % get gt image file name by k
    GT = imread(fullfile(pwd, sprintf('IMG_%02d_GT.png', k)));
    targets = GT(:,:,1);
    if flag
        % If you only have background and bacteria, you need to change the flag to 0, 1
        targets(targets==2) = 1;
    end
    % get precision,recall
    [~,precision,recall] = bfscore(double(outputs),double(targets));
    % get dice score
    score = dice(double(outputs),double(targets));
    score_matrix_k = [score,precision,recall];
    score_matrix_k = round(score_matrix_k*10000)/10000;
    score_matrix1 = score_matrix_k(1,:);
    score_matrix2 = [];
    if size(score_matrix_k, 1) > 1
        % If there is background, cells and bacteria, the evaluation results for bacteria are included
        score_matrix2 = score_matrix_k(2,:);
    end
    s = [];
    % get dice_score, precision, recall
    s(1).name = "cells";
    s(1).dice_score = score_matrix1(1);
    s(1).precision = score_matrix1(2);
    s(1).recall = score_matrix1(3);
    s(2).name = "bacteria";
    if ~isempty(score_matrix2)
        % If there is background, cells and bacteria, the evaluation results for bacteria are included
        s(2).dice_score = score_matrix2(1);
        s(2).precision = score_matrix2(2);
        s(2).recall = score_matrix2(3);
    end% all resuilt will be print in command windows ,pleas check it 
    T = struct2table(s)
    
    % To visualise the ground truth image, you can
    % use the following code.
    L_GT = label2rgb(GT(:,:,1), 'prism','k','shuffle');
    figure, imshow(L_GT); title('GT image');
    
    close all
    for i = 1 : length(s)
        data{end+1,1} = k;
        data{end,2} = s(i).name;
        data{end,3} = s(i).dice_score;
        data{end,4} = s(i).precision;
        data{end,5} = s(i).recall;
    end
end
% pirnt all result of 15 images in file named "data.xlsx"
xlswrite('data.xlsx',data)


%-----------------------------------------------------------------------
%To improve the aesthetics of the code
%To improve code readability
%To reduce code re-use
%I encapsulate the processing of segment in sub-functions
% i hvae test put all subfuntion in task 5 but faill
%sub function must be placed after all functions
%--------------------------------------------------------------

%% sub functions
%---------------------------------------------------------------------------
% segment cells main function
function res = main_seg1(I)
% first seg,use rate 1, normal rate :Processing bright red blood cells
bw1 = get_bw1(I,1);
% second seg,use rate 0.25, normal rate:Processing dark red blood cells 
bw2 = get_bw1(I,0.25);
% clear by first ,In dark cells, remove bright cells and keep only dark cells
[r,c] = find(bw1);
bw3 = bwselect(bw2, c, r);
bw2(bw3) = 0;
% filter noise Remove noise from dark cells
bw2 = filter_bw1(bw2);
% get result: Merge light and dark cells to get a complete image of the cell
res = logical(bw1+bw2);
end

%---------------------------------------------------------------------------
% segment bacteria  main function
function res2 = main_seg2(I,res)
% first seg,use rate 1, normal rate deal with bright bacteria
bw1 = get_bw2(I,res,1);
% second seg,use rate 0.3 Deal with dark bacteria
bw2 = get_bw2(I,res,0.3);
% clear by first£¬In dark bacterial images, eliminate areas of light
% bacteria
[r,c] = find(bw1);
bw3 = bwselect(bw2, c, r);
bw2(bw3) = 0;
[~,num] = bwlabel(bw1);
if num == 3
    % Eliminate the noise of dark bacteria
    bw1 = filter_bw2(bw1);
end
% get result
bw2 = filter_bw2(bw2);
% Merge light and dark bacteria into one complete image
res2 = logical(bw1+bw2);
[r,c] = find(res);
bw3 = bwselect(res2, c, r);
if ~isequal(bw3,res2)
    % Eliminate red blood cell interference
    res2(bw3) = 0;
    if length(find(res2))/length(find(bw3)) < 0.3
        res2 = bw3;
    end
end
end

%---------------------------------------------------------------------------
% filter cell noise
function openBW = filter_bw1(openBW)
% open to clear small area object
openBW = bwareaopen(openBW,300);
[L,num] = bwlabel(openBW);
% label props
stats = regionprops(L,'Area','BoundingBox','Solidity');
for i = 1 : num
    % delete noise object
    %Keep targets with large areas and high infill
    if stats(i).Area > 4000 && stats(i).Solidity > 0.9
    else
        openBW(L==i) = 0;
    end
end
% delete line object
openBW = imopen(openBW, strel('disk', 9));
end


%---------------------------------------------------------------------------
% filter line noise
function openBW = filter_bw2(openBW)
% clear border noise
openBW = imclearborder(openBW);
% clear small area
openBW = bwareaopen(openBW,300);
[L,num] = bwlabel(openBW);
stats = regionprops(L,'Area','BoundingBox','Solidity');
for i = 1 : num
    % Filter by area, retain blobs whose area and solidity meet the range, and delete others
    if (((stats(i).Area > 100 && stats(i).Area < 700)||(stats(i).Area > 1500))...
            && stats(i).Solidity > 0.5 && stats(i).Solidity < 0.7) || ...
            (stats(i).Area > 700 && stats(i).Solidity > 0.9 && max(stats(i).BoundingBox(3:4)) < 80)
    else
        openBW(L==i) = 0;
    end
end
end

%---------------------------------------------------------------------------
% segment cell by thresholds 
% Use different thresholds to segment cells, with normal thresholds for bright cells and reduced thresholds for dark cells
function openBW = get_bw1(I,rate)
if nargin < 2
    % thresh rate for segment
    rate = 1;
end
% 1.grayscale image
I_grayscale = rgb2gray(I);
%2.resised image
I_rescale = imresize(I_grayscale,[512 NaN]);
%3. enhanced image
J = 255*im2double(I_rescale);
mi = min(min(J));
ma = max(max(J));
I_Enhance = imadjust(I_rescale,[mi/255; ma/255],[0 1]);
%4. binary image
thresh=graythresh(I_Enhance);
% use rate to threshold
I_bina = imbinarize(I_Enhance,thresh*rate);
if thresh > 0.4 && rate == 1
    % use for different process
    I_bina = imbinarize(I_Enhance,thresh*0.8);
    I_bina = imfill(I_bina,'holes');
end
%5.Edge detection
I_edge = edge(I_bina,'canny',0.99);
SE= strel('disk',1);
BWsdil = imdilate(I_edge,SE);
%Some cells are close to the edge, and the calculation area will be affected, 
% to prevent this, draw lines on the top, bottom, left, and right edges of the image in advance.
col1 = find(BWsdil(1, :), 1, 'first');
col2 = find(BWsdil(1, :), 1, 'last');
BWsdil(1, col1:col2) = true;
col1 = find(BWsdil(end, :), 1, 'first');
col2 = find(BWsdil(end, :), 1, 'last');
BWsdil(end, col1:col2) = true;
row1 = find(BWsdil(:, end), 1, 'first');
row2 = find(BWsdil(:, end), 1, 'last');
BWsdil(row1:row2, end) = true;
col1 = find(BWsdil(:, 1), 1, 'first');
col2 = find(BWsdil(:, 1), 1, 'last');
BWsdil( col1:col2,1 ) = true;
BWsdil = imfill(BWsdil, 'holes');
if thresh >= 0.42 && rate == 1
    BWsdil=logical(BWsdil+I_bina);
    % denoise small block
    BWsdil=bwareaopen(BWsdil,300);
    [~,num] = bwlabel(BWsdil);
    if num <= 5
        % filter block link
        BWsdil=imerode(BWsdil,strel('disk', num));
    else
        BWsdil=imerode(BWsdil,strel('disk', 9));
    end
end
SE2 = strel('disk',1);
openBW = imclose(BWsdil,SE2);
openBW = bwareaopen(openBW,300);
[L,num] = bwlabel(openBW);
stats = regionprops(L,'Area','BoundingBox','Solidity');
for i = 1 : num
    % use for region property
    %Calculate area and solidity and aspect ratio, remove interference areas
    if stats(i).Area < 3000 && stats(i).Solidity < 0.9 && max(stats(i).BoundingBox(3:4))/min(stats(i).BoundingBox(3:4)) > 1.5
        openBW(L==i) = 0;
    end%Determine the area and solidity,
    % calculate the maximum length and width, and delete the interference area
    if stats(i).Area < 500 && max(stats(i).BoundingBox(3:4)) < 50
        openBW(L==i) = 0;
    end
end
[~,num] = bwlabel(openBW);
if num <= 10
    % use for multi block process
    openBW = imopen(openBW, strel('disk', 9));
    [L,num] = bwlabel(openBW);
    % get stat property
    stats = regionprops(L,'Area','BoundingBox','Solidity');
    for i = 1 : num
        if stats(i).Area < 500 && max(stats(i).BoundingBox(3:4)) < 50
            %Remove distractions again based on area and BoundingBox
            openBW(L==i) = 0;
        end
    end
else
    % use for area filter
    I_bina2 = imbinarize(I_Enhance,thresh*rate/2);
    I_bina2 = imfill(I_bina2,'holes');%fill all remaining areas
    b1 = bwareafilt(openBW, 1);
    b2 = bwareafilt(openBW, 1, 'smallest');
    % process big and small area
    [r,c] = find(b1); c1 = bwselect(I_bina2,c,r);
    [r,c] = find(b2); c2 = bwselect(I_bina2,c,r);
    c2 = imclose(c2, strel('disk', 100));
    if numel(find(c1))/numel(find(I_bina2)) > 0.3
        % can be add to image
        openBW = logical(openBW+c1+c2);
    end
end
[~,num] = bwlabel(openBW);
if num < 3
    % use for denoise
    openBW = bwareaopen(openBW,3000);
end
end


%---------------------------------------------------------------------------
% segment bacteria: Segment bacteria using different thresholds
function bt = get_bw2(I,res,rate)
if nargin < 3
    rate = 1;
end
% 1.grayscale image
I_grayscale = rgb2gray(I);%cover image to graysacle
%figure,imshow(I_grayscale);title('grayscale image');
%2.resised image
I_rescale = imresize(I_grayscale,[512 NaN]);
%figure,imshow(I_rescale);title('resized Image ');
%3. enhanced image
J = 255*im2double(I_rescale);
mi = min(min(J)); %get minimum value
ma = max(max(J));% get maximunm value
I_Enhance = imadjust(I_rescale,[mi/255; ma/255],[0 1]);
%4. binary image
thresh=graythresh(I_Enhance);
if rate == 1
    % reduce cell normal
    I_Enhance(res) = 0;
    I_bina = imbinarize(I_Enhance,thresh);
    I_bina(res) = 0;
else
    % reduce cell by dilate
    res = imdilate(res, strel('disk', 15));
    I_Enhance(res) = 0;
    I_bina = imbinarize(I_Enhance,thresh*rate);
    I_bina(res) = 0;
end
%5.Edge detection
I_edge = edge(I_bina,'canny',0.99);
% 6.Image segmentation
SE= strel('disk',1);
BWsdil = imdilate(I_edge,SE);
% Make top line white:
col1 = find(BWsdil(1, :), 1, 'first');
col2 = find(BWsdil(1, :), 1, 'last');
BWsdil(1, col1:col2) = true;
% Make bottom line white:
col1 = find(BWsdil(end, :), 1, 'first');
col2 = find(BWsdil(end, :), 1, 'last');
BWsdil(end, col1:col2) = true;
% Make right line white:
row1 = find(BWsdil(:, end), 1, 'first');
row2 = find(BWsdil(:, end), 1, 'last');
BWsdil(row1:row2, end) = true;
% Make left line white:
col1 = find(BWsdil(:, 1), 1, 'first');
col2 = find(BWsdil(:, 1), 1, 'last');
BWsdil( col1:col2,1 ) = true;
BWsdil = imfill(BWsdil, 'holes');
SE2 = strel('disk',1);% Create a circle with radius 1 as the structural element
openBW = imclose(BWsdil,SE2);
[L,num] = bwlabel(openBW);
stats = regionprops(L);
for i = 1 : num
    % filter object area
    if stats(i).Area < 200 && max(stats(i).BoundingBox(3:4)) < 10
        % Remove interference areas based on area, solidity and maximum length and width
        openBW(L==i) = 0;
    end
end
[B,L] = bwboundaries(openBW,4);
L2 = L;
stats = regionprops(L,'Area','Centroid');
props = regionprops(L, 'MajorAxisLength', 'MinorAxisLength');
threshold = 0.7;
bt = logical(zeros(size(I_bina)));
for k = 1:length(B)
    boundary = B{k};% get coordinates of label k
    delta_sq = diff(boundary).^2; %calculate perimeter
    perimeter = sum(sqrt(sum(delta_sq,2)));
    area = stats(k).Area;% get area of label k
    metric = 4*pi*area/perimeter^2;%Compute the roundness metric
    metric_string = sprintf('%.2f',metric);%print result
    aspectRatios = [props(k).MajorAxisLength] ./ [props(k).MinorAxisLength];
    conditon = aspectRatios >=2;
    if conditon == 1
        % use for line
        if metric < threshold && metric > 0.05
            bt(L2 == k) = 1;
            bt = logical(bt);
        end
    end
end
bt = logical(bt);
end


