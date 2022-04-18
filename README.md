## Task 1 --Pre-processing
Load the image into the matlab script using the function imread(),Convert the image to grayscale using the function rgb2gray().Use the function resize() to modify the size of the image, the height is 512, NaN is used as the width, and matlab will automatically adjust to the same aspect ratio. Use the imhist() function to draw a histogram.

Discuss: Image enhancement techniques include contrast stretching and histogram equalisation. The imadjust() function is used to extend the contrast. The histq() and adapthisteq() functions are used in histogram equalisation. After testing with many images, the results returned by histq() and adapthisteq() enhance the noise, indicating that histogram equalisation is not appropriate. Use the imadjust() function and the contrast stretch technique. Set the gamma value to 0.7 after some try; a gamma less than 0.6 will accentuate the noise. The original visual detail of cricle edge will vanish above 1.0. The rest of the settings are set the same way they were in the workshop. The data that was before concentrated in 25 is now dispersed over the range of 0 to 200 as a result of the expanded histogram observation.

## Task2–Edge Detection
In MATLAB, the edge() function is used to detect the edge. The edge() function has a total of seven methods. The Canny technique was adopted after several experiments. Canny's method involves determining the local maximum of object's gradient. The gradient is computed using the Gaussian filter's derivative by the edge function. This method detects both strong and weak edges using a double threshold and includes the weak edge in the output if it is connected to the strong edge. Canny's approach is less susceptible to noise than other methods and is more likely to identify genuinely weak edges due to the use of dual thresholds. The Canny method has three advantages: 1. The Canny operator is unaffected by noise. 2. The existence of a Gaussian filter enables the picture to be free of noise; 3. Using a threshold approach, edges may be detected in a noisy condition.

## Task3– Simple Segmentation
Due to the proximity of the red blood cells to the picture margins, I opted out of employing morphological processing approaches (imdilate, imerode, imopen and imclose). Select the activecontour() function, which was taught in lecture week 9. To begin, create an initial contour (mask) around the item and segment the picture using the activecontour() function. By default, the function iteratively develops the split. However, the picture remains unusable, and the number of iterations is raised to 500. Finally, the picture is divided successfully.

## Task4–Object Recognition
Solution: 1. Trace the exterior bounds of objects using the bwboundaries() method. Make a note of the outer boundaries of each blob and use the plot() function to print out all of the blobs' boundaries. 2. Using the regionprops() function, determine the area and centroid of the image's object. Second, a threshold value of 0.33 is used as the criterion for determining if the figure is circular. Additionally, get the outside border of each item included inside the loop and compute its perimeter, area, and roundness. Finally, the image is evaluated circular according to the previously established threshold, and the label for the image identified as a circle is set to 2, while the label for the non-circular image is put to 3. 3. Convert the label matrix to an RCG image using the label2rgb() function, and change the backdrop colour to black

## Task5–Robust Method
In this algorithm, there are 3 steps to do: 1. Preprocessing. 2. Segment the image. 3. Mark and color images. I chose IMG_04 as the example image.

### 1.Preprocessing:
There are three steps in preprocessing. 1. Use the imread() function to read the file. 2. Use the rgb2gray() function to convert the read file to a grayscale image. 3. Use the imadjust() function to enhance the image. In the imadjust() function, the input threshold is set to [mi/255; ma/255], and the output threshold is set to [0 1]. (mi = min(min(The original image)), ma = max (max (The original image))). It is benefit to enhance image that get more detail about bacteria and red cell.

### 2.Segment the image:
To cleanly divide the image. I performed two divides in the "fragment Image" phase. Using the method of reorganizing the image after segment the image twice can effectively remove noise and avoid setting too large a mask area to shield small red blood cells and bacteria. This method helps to improve the final accuracy
The first step is to segment the figure's red blood cells. Red blood cells are separated into bright and dark cells in each of the 15 images. For bright red blood cells, a normal threshold (threshold of 1) was employed; for dark red blood cells, a decreased threshold (threshold of 0.25) was utilised. When splitting dark red blood cells, the bright red blood cells must be excluded, which is accomplished using the bwselect() function. After obtaining the dark red blood cells, noise reduction is conducted, and the processed light red blood cells and dark red blood cells are combined into a single picture using the logical() function to create a full image. The process of splitting red blood cells and removing noise is enclosed in sub-functions; due to a lack of space for a detailed explanation, please refer to the code.

Segmenting the bacteria in the image is the second step. Bacteria are split into bright and dark colours in the 15 photos. For light-colored bacteria, use a standard threshold (threshold of 1) and for dark-colored bacteria, use a decreased threshold (threshold of 0.3). Following the removal of bright bacteria using the bwselect() function, both dark and bright bacteria are segment and noise is eliminated independently. Finally, combine the processed brilliant red blood cells and dark red blood cells into a single picture using the logical() method to get a full image. The process of splitting red blood cells and removing noise is enclosed in sub-functions，please refer to the code.

### 3.Mark and color images:
Label red blood cells, bacteria and background as 1, 2, and 0, respectively, as required by the assessment. Because the background is 0, it is directly ignored. Convert the image to an RBG image using the label2rgb() function.

### Task6–Performance Evaluation
The dice score is calculated using the dice () function, the precision and recall are calculated using the bfscore() function .

### Data description:
According to the statistics, dice score, precision, and recall all have extremely high values, indicating that the object recognition impact is rather excellent. The standard deviations are all quite small, implying that the variations are also quite minor. All the results are recorded in the file named "data.xlsx", you can see the specific score of each image.










