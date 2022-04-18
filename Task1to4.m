clear; close all;

% Task 1: Pre-processing -----------------------
% Step-1: Load input image
I = imread('IMG_11.png');
%figure, imshow(I)

% Step-2: Covert image to grayscale
I_gray = rgb2gray(I);
%figure, imshow(I_gray)

% Step-3: Rescale image
Rescale_img = imresize(I_gray,[512 NaN]);%set image high is 512,use Nan can change image depends on aspect ratio
%figure,imshow(Rescale_img);title('resize image');

% Step-4: Produce histogram before enhancing
%figure,imhist(Rescale_img,64);

% Step-5: Enhance image before binarisation
J = 255*im2double(Rescale_img);
mi = min(min(J)); %get minimum value
ma = max(max(J));% get maximunm value
Enhance_img = imadjust(Rescale_img,[mi/255; ma/255],[0 1]);%test differ value in function
figure,imshow(Enhance_img);title('Enhanced image');

% Step-6: Histogram after enhancement
%figure,imhist(Enhance_img,64);

% Step-7: Image Binarisation
thresh=graythresh(Enhance_img);%get threshold value
Binarisation_img = imbinarize(Enhance_img,thresh);
%figure,imshow(Bina_img);title('Binary figure');

% Task 2: Edge detection ------------------------
Edge_img = edge(Binarisation_img,'canny',0.86);%select Canny as method to get edge
figure,imshow(Edge_img);title('edge detection of Canny');

% Task 3: Simple segmentation --------------------
mask = zeros(size(Enhance_img));%Specifies the initial outline around the object of interest
mask(25:end-25,25:end-25) = 1;
BW = activecontour(Enhance_img,mask,500);%Segment the image using the activecontour function. By default, the function evolves the split over 500 iterations
closeBW = bwareaopen(BW, 700);%Remove white blobs with area less than 700
Fill_img = imfill(closeBW,'holes');% fill all holes
figure,imshow(Fill_img);title('segmentation');


% Task 4: Object Recognition --------------------
[B,L] = bwboundaries(Fill_img,4);% get number of bwboundaries
L_copy = L;
%figure,imshow(label2rgb(L,@jet,[.0 .0 .0]));
hold on

stats = regionprops(L,'Area','Centroid');% get area and centroid of blobs
% set a threshold is 0.33
threshold = 0.33;

for k = 1:length(B)
    boundary = B{k};% get coordinates of label k
    delta_sq = diff(boundary).^2; %calculate perimeter 
    perimeter = sum(sqrt(sum(delta_sq,2)));
    area = stats(k).Area;% get area of label k
    metric = 4*pi*area/perimeter^2;%Compute the roundness metric
    metric_string = sprintf('%2.2f',metric);%print result

    if metric > threshold % According to Shape Factor, judge cell or bacteria
    L(L_copy == k) = 3;% label cells
    else
    L(L_copy == k) = 2;% label bactria 
    end
% Add a text description of the result
    text(boundary(1,2)-50,boundary(1,1)+13,metric_string,'Color','r','FontSize',14,'FontWeight','bold')
end
%show result
result = label2rgb(L,@jet,[.0 .0 .0]);
figure,imshow(result);


