// Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include "math.h"
#include <queue>
#include <iostream>
#include <utility>
using namespace std;


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
: pixels(NULL),
npixels(0),
width(0),
height(0)
{
}



R2Image::
R2Image(const char *filename)
: pixels(NULL),
npixels(0),
width(0),
height(0)
{
	// Read image
	Read(filename);
}



R2Image::
R2Image(int width, int height)
: pixels(NULL),
npixels(width * height),
width(width),
height(height)
{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
: pixels(NULL),
npixels(width * height),
width(width),
height(height)
{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++)
		pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
: pixels(NULL),
npixels(image.npixels),
width(image.width),
height(image.height)

{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++)
		pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
	// Free image pixels
	if (pixels) delete[] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
	// Delete previous pixels
	if (pixels) { delete[] pixels; pixels = NULL; }

	// Reset width and height
	npixels = image.npixels;
	width = image.width;
	height = image.height;

	// Allocate new pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++)
		pixels[i] = image.pixels[i];

	// Return image
	return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2, 3.5);
	R2Point p2(2.1, 2.2);
	R2Point p3(0.2, 1.6);
	R2Point p4(0.0, 0.5);
	R2Point p5(-0.2, 4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1, 5, 1, 6);

	linEquations[1][1] = p1[0] * p1[0];
	linEquations[1][2] = p1[0] * p1[1];
	linEquations[1][3] = p1[1] * p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0] * p2[0];
	linEquations[2][2] = p2[0] * p2[1];
	linEquations[2][3] = p2[1] * p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0] * p3[0];
	linEquations[3][2] = p3[0] * p3[1];
	linEquations[3][3] = p3[1] * p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;

	linEquations[4][1] = p4[0] * p4[0];
	linEquations[4][2] = p4[0] * p4[1];
	linEquations[4][3] = p4[1] * p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0] * p5[0];
	linEquations[5][2] = p5[0] * p5[1];
	linEquations[5][3] = p5[1] * p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n", p1[0], p1[1]);
	printf("Point #2: %f,%f\n", p2[0], p2[1]);
	printf("Point #3: %f,%f\n", p3[0], p3[1]);
	printf("Point #4: %f,%f\n", p4[0], p4[1]);
	printf("Point #5: %f,%f\n", p5[0], p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1, 6, 1, 6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n", singularValues[1], singularValues[2], singularValues[3], singularValues[4], singularValues[5], singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for (int i = 2; i<7; i++) if (singularValues[i]<singularValues[smallestIndex]) smallestIndex = i;
	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n", nullspaceMatrix[1][smallestIndex], nullspaceMatrix[2][smallestIndex], nullspaceMatrix[3][smallestIndex], nullspaceMatrix[4][smallestIndex], nullspaceMatrix[5][smallestIndex], nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n", p1[0] * p1[0] * nullspaceMatrix[1][smallestIndex] +
		p1[0] * p1[1] * nullspaceMatrix[2][smallestIndex] +
		p1[1] * p1[1] * nullspaceMatrix[3][smallestIndex] +
		p1[0] * nullspaceMatrix[4][smallestIndex] +
		p1[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n", p2[0] * p2[0] * nullspaceMatrix[1][smallestIndex] +
		p2[0] * p2[1] * nullspaceMatrix[2][smallestIndex] +
		p2[1] * p2[1] * nullspaceMatrix[3][smallestIndex] +
		p2[0] * nullspaceMatrix[4][smallestIndex] +
		p2[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n", p3[0] * p3[0] * nullspaceMatrix[1][smallestIndex] +
		p3[0] * p3[1] * nullspaceMatrix[2][smallestIndex] +
		p3[1] * p3[1] * nullspaceMatrix[3][smallestIndex] +
		p3[0] * nullspaceMatrix[4][smallestIndex] +
		p3[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n", p4[0] * p4[0] * nullspaceMatrix[1][smallestIndex] +
		p4[0] * p4[1] * nullspaceMatrix[2][smallestIndex] +
		p4[1] * p4[1] * nullspaceMatrix[3][smallestIndex] +
		p4[0] * nullspaceMatrix[4][smallestIndex] +
		p4[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n", p5[0] * p5[0] * nullspaceMatrix[1][smallestIndex] +
		p5[0] * p5[1] * nullspaceMatrix[2][smallestIndex] +
		p5[1] * p5[1] * nullspaceMatrix[3][smallestIndex] +
		p5[0] * nullspaceMatrix[4][smallestIndex] +
		p5[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34, -2.8);

	printf("A point off the conic: %f\n", test_point[0] * test_point[0] * nullspaceMatrix[1][smallestIndex] +
		test_point[0] * test_point[1] * nullspaceMatrix[2][smallestIndex] +
		test_point[1] * test_point[1] * nullspaceMatrix[3][smallestIndex] +
		test_point[0] * nullspaceMatrix[4][smallestIndex] +
		test_point[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
	// Brighten the image by multiplying each pixel component by the factor.
	// This is implemented for you as an example of how to access and set pixels
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) *= factor;
			Pixel(i, j).Clamp();
		}
	}
}

void R2Image::
SobelX(void)
{
	// Apply the Sobel oprator to the image in X direction
	R2Image temp(width, height);

	// Sobel Kernel
	float sobel_x[3][3] = {
		{ 1, 0, -1 },
		{ 2, 0, -2 },
		{ 1, 0, -1 } };

	// Applying the mask to every pixel in the image and saves it on a new temp image
	for (int i = 1; i <= width - 2; i++) {
		for (int j = 1; j <= height - 2; j++) {
			R2Pixel pix_x = (sobel_x[0][0] * Pixel(i - 1, j - 1)) + (sobel_x[0][1] * Pixel(i, j - 1)) + (sobel_x[0][2] * Pixel(i + 1, j - 1)) +
				(sobel_x[1][0] * Pixel(i - 1, j)) + (sobel_x[1][1] * Pixel(i, j)) + (sobel_x[1][2] * Pixel(i + 1, j)) +
				(sobel_x[2][0] * Pixel(i - 1, j + 1)) + (sobel_x[2][1] * Pixel(i, j + 1)) + (sobel_x[2][2] * Pixel(i + 1, j + 1));
			temp.Pixel(i, j) = pix_x;
			//temp.Pixel(i, j).Clamp();
		}
	}

	// Copies the values of the temp image pixels onto the original image
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) = temp.Pixel(i, j);
		}
	}
}

void R2Image::
SobelY(void)
{
	// Apply the Sobel oprator to the image in Y direction
	R2Image temp(width, height);
	float sobel_y[3][3] = {
		{ -1, -2, -1 },
		{ 0, 0, 0 },
		{ 1, 2, 1 } };

	// Applying the mask to every pixel in the image and saves it on a new temp image
	for (int i = 1; i <= width - 2; i++) {
		for (int j = 1; j <= height - 2; j++) {
			R2Pixel pix_y = (sobel_y[0][0] * Pixel(i - 1, j - 1)) + (sobel_y[0][1] * Pixel(i, j - 1)) + (sobel_y[0][2] * Pixel(i + 1, j - 1)) +
				(sobel_y[1][0] * Pixel(i - 1, j)) + (sobel_y[1][1] * Pixel(i, j)) + (sobel_y[1][2] * Pixel(i + 1, j)) +
				(sobel_y[2][0] * Pixel(i - 1, j + 1)) + (sobel_y[2][1] * Pixel(i, j + 1)) + (sobel_y[2][2] * Pixel(i + 1, j + 1));
			temp.Pixel(i, j) = pix_y;
			//temp.Pixel(i, j).Clamp();
		}
	}

	// Copies the values of the temp image pixels onto the original image
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) = temp.Pixel(i, j);
		}
	}
}

void R2Image::
LoG(void)
{
	// Apply the LoG oprator to the image

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
	// Changes the saturation of an image
	// Find a formula that changes the saturation without affecting the image brightness

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////

void R2Image::
Sharpen()
{
	// Sharpen Kernel
	int weight[3][3] = { 
		{ 1, 1, 1 },
		{ 1, -8, 1 },
		{ 1, 1, 1 } };

	double factor = 0.2;

	R2Image copy = R2Image(*this);
	R2Pixel pix, newPix;
	int newVal;
	R2Pixel nullPix = R2Pixel(0.0, 0.0, 0.0, 1.0);

	for (int dy = 1; dy < height - 1; dy++) {
		for (int dx = 1; dx < width - 1; dx++) {
			pix = R2Pixel();
			for (int j = -1; j < 2; j++) {
				for (int i = -1; i < 2; i++) {
					pix += weight[j + 1][i + 1] * Pixel(dy + j, dx + i);
				}
			}
			newPix = Pixel(dy, dx) - factor * pix;
			newVal = newPix.Red() + newPix.Green() + newPix.Blue();
			if (newVal < 0) newVal = 0;
			Pixel(dy, dx) = newPix;
		}
	}
}


void R2Image::
Blur(double sigma)
{
	R2Image temp(*this);

	double FACTOR = 1.0 / sqrt(2 * 3.14*sigma*sigma); //constant dependent on sigma
	int kernelSize = (int)(6 * sigma) + 1;
	float* kernel = (float*)malloc(kernelSize * sizeof(float));
	double totalWeight = 0;

	// build kernel
	int w = (int)(kernelSize - 1) / 2;
	for (int i = -1 * w; i <= w; i++){
		kernel[i + w] = FACTOR * exp(((-1 * i*i) / (2 * sigma*sigma)));
		totalWeight += kernel[i + w];

	}

	// normalize kernel
	for (int i = 0; i<kernelSize; i++){
		kernel[i] = kernel[i] / totalWeight;
	}

	// temporary pixel to be used in intermediate computation
	R2Pixel val;

	// vertical pass
	for (int x = 0; x < width; x++){
		for (int y = w; y < height - w; y++){

			// compute pixel value
			val = R2Pixel();
			for (int k = -w; k <= w; k++){
				val += (kernel[k + w])*Pixel(x, y + k);
			}
			// set pixel
			temp.SetPixel(x, y, val);
		}
	}

	// horizontal pass
	for (int y = 0; y < height; y++){
		for (int x = w; x < width - w; x++){

			// compute pixel value
			val = R2Pixel();
			for (int k = -w; k <= w; k++){
				val += (kernel[k + w])*temp.Pixel(x + k, y);
			}
			// set pixel
			SetPixel(x, y, val);
		}
	}

}


struct Pix {
	R2Pixel pix;
	int x;
	int y;
} typedef pixCompare;

struct Compare {
	bool  operator() (pixCompare &left, pixCompare &right)
	{
		return left.pix.Luminance() < right.pix.Luminance();
	}
};



priority_queue<pixCompare, vector<pixCompare>, Compare> MakePriorityQueue(R2Image &img)
{
	priority_queue<pixCompare, vector<pixCompare>, Compare> pq;


	int i, j;
	for (i = 0; i < img.Width(); i++) {
		for (j = 0; j < img.Height(); j++) {
			pixCompare addedPix;
			addedPix.pix = img.Pixel(i, j);
			addedPix.x = i;
			addedPix.y = j;
			pq.push(addedPix);
			img.Pixel(i, j).Clamp();
		}
	}


	return pq;
}



void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
	if (x0>x1)
	{
		int x = y1;
		y1 = y0;
		y0 = x;

		x = x1;
		x1 = x0;
		x0 = x;
	}
	int deltax = x1 - x0;
	int deltay = y1 - y0;
	float error = 0;
	float deltaerr = 0.0;
	if (deltax != 0) deltaerr = fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
	// note that this division needs to be done in a way that preserves the fractional part
	int y = y0;
	for (int x = x0; x <= x1; x++)
	{
		Pixel(x, y).Reset(r, g, b, 1.0);
		error = error + deltaerr;
		if (error >= 0.5)
		{
			if (deltay>0) y = y + 1;
			else y = y - 1;

			error = error - 1.0;
		}
	}
	if (x0>3 && x0<width - 3 && y0>3 && y0<height - 3)
	{
		for (int x = x0 - 3; x <= x0 + 3; x++)
		{
			for (int y = y0 - 3; y <= y0 + 3; y++)
			{
				Pixel(x, y).Reset(r, g, b, 1.0);
			}
		}
	}
}


void FeatureDetectionOld(R2Image &img, priority_queue<pixCompare, vector<pixCompare>, Compare> pq) {

	int COUNT = 150;
	int selectedSoFar = 0;
	int pixelsApart = 10;
	R2Image blank = R2Image(img.Width(), img.Height());
	R2Pixel whitePix = R2Pixel(1.0, 1.0, 1.0, 1);
	R2Pixel blackPix = R2Pixel(0.0, 0.0, 0.0, 1);
	R2Pixel redPix = R2Pixel(1.0, 0.0, 0.0, 1);

	int i, j;
	for (i = 0; i < img.Width(); i++) {
		for (j = 0; j < img.Height(); j++) {
			blank.Pixel(i, j) = whitePix;
		}
	}
	pixCompare currPix;
	int x_coord, y_coord;

	double radius = 1.0;

	while (selectedSoFar < COUNT) {
		currPix = pq.top();
		pq.pop();
		x_coord = currPix.x;
		y_coord = currPix.y;
		if (!blank.Pixel(x_coord, y_coord).IsBlack()) {
			selectedSoFar++;
			for (int dx = max(0, x_coord - pixelsApart); dx < min(blank.Width(), x_coord + pixelsApart); dx++) {
				for (int dy = max(0, y_coord - pixelsApart); dy < min(blank.Height(), y_coord + pixelsApart); dy++) {
					blank.Pixel(dx, dy) = blackPix;
				}
			}

			for (int dx = max(0, x_coord - pixelsApart); dx < min(img.Width(), x_coord + pixelsApart); dx++) {
				printf("dx = %d \n", dx);
				img.Pixel(dx, max(0, y_coord - pixelsApart)) = redPix;
				img.Pixel(dx, min(img.Height(), y_coord + pixelsApart)) = redPix;
			}
			for (j = max(0, y_coord - pixelsApart); j < min(img.Height(), y_coord + pixelsApart); j++) {
				img.Pixel(max(0, x_coord - pixelsApart), j) = redPix;
				img.Pixel(min(img.Width(), x_coord + pixelsApart), j) = redPix;
			}

			radius += 0.06;
		}
	}
}


pixCompare* FeatureDetectionNew(R2Image &img, R2Image &original, priority_queue<pixCompare, vector<pixCompare>, Compare> pq) {

	int COUNT = 150;
	int selectedSoFar = 0;
	int pixelsApart = 10;

	R2Image blank = R2Image(img.Width(), img.Height());
	R2Pixel whitePix = R2Pixel(1.0, 1.0, 1.0, 1);
	R2Pixel blackPix = R2Pixel(0.0, 0.0, 0.0, 1);
	//R2Pixel redPix = R2Pixel(1.0, 0.0, 0.0, 1);

	pixCompare* featuresArray = new pixCompare[COUNT];
	int index = 0;


	int i, j;
	for (i = 0; i < blank.Width(); i++) {
		for (j = 0; j < blank.Height(); j++) {
			blank.Pixel(i, j) = whitePix;
		}
	}

	pixCompare currPix;
	int x_coord, y_coord;


	while (selectedSoFar < COUNT) {
		currPix = pq.top();
		pq.pop();
		x_coord = currPix.x;
		y_coord = currPix.y;
		if (!blank.Pixel(x_coord, y_coord).IsBlack()) {
			featuresArray[index] = currPix;
			index++;
			selectedSoFar++;
			for (int dx = max(0, x_coord - pixelsApart); dx < min(blank.Width(), x_coord + pixelsApart); dx++) {
				for (int dy = max(0, y_coord - pixelsApart); dy < min(blank.Height(), y_coord + pixelsApart); dy++) {
					blank.Pixel(dx, dy) = blackPix;
				}
			}
		}
	}

	return featuresArray;
}


void R2Image::
Harris(double sigma)
{
	// Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

	R2Image copy = R2Image(*this);

	R2Pixel offset(0.5, 0.5, 0.5, 1);
	R2Image sobel_X(*this);
	sobel_X.SobelX();
	R2Image sobel_Y(*this);
	sobel_Y.SobelY();
	R2Image sobel_XY(*this);


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			sobel_XY.Pixel(i, j) = sobel_X.Pixel(i, j) * sobel_Y.Pixel(i, j);
			sobel_X.Pixel(i, j) = sobel_X.Pixel(i, j) * sobel_X.Pixel(i, j);
			sobel_Y.Pixel(i, j) = sobel_Y.Pixel(i, j) * sobel_Y.Pixel(i, j);
		}
	}

	sobel_X.Blur(sigma);
	sobel_Y.Blur(sigma);
	sobel_XY.Blur(sigma);

	R2Pixel tempPix;
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			tempPix = sobel_X.Pixel(x, y) * sobel_Y.Pixel(x, y) - sobel_XY.Pixel(x, y) * sobel_XY.Pixel(x, y) -
				0.04 * ((sobel_X.Pixel(x, y) + sobel_Y.Pixel(x, y)) * (sobel_X.Pixel(x, y) + sobel_Y.Pixel(x, y))) + offset;

			copy.Pixel(x, y) = tempPix;
		}
	}

	
	priority_queue<pixCompare, vector<pixCompare>, Compare> pq = MakePriorityQueue(copy);

	// get the coordinates of the top features
	FeatureDetectionOld(*this, pq);
}




double singleSSD(int x, int y, int x_start, int y_start, R2Image &img, R2Image &other) {

	double result = 0.0;

	double diffRed = 0.0;
	double diffGreen = 0.0;
	double diffBlue = 0.0;

	int featureSpan = 10;

	int i, j;

	for (i = -1 * featureSpan; i <= featureSpan; i++) {
		for (j = -1 * featureSpan; j <= featureSpan; j++) {
			if (x_start + i >= 0 && x_start + i < img.Width() &&
				y_start + j >= 0 && y_start + j < img.Height() &&
				x + i >= 0 && x + i < img.Width() &&
				y + j >= 0 && y + j < img.Height()) {


				diffRed = img.Pixel(i + x, j + y).Red() - other.Pixel(i + x_start, j + y_start).Red();
				diffGreen = img.Pixel(i + x, j + y).Green() - other.Pixel(i + x_start, j + y_start).Green();
				diffBlue = img.Pixel(i + x, j + y).Blue() - other.Pixel(i + x_start, j + y_start).Blue();

				result += (diffRed * diffRed) + (diffGreen * diffGreen) + (diffBlue * diffBlue);
			}
		}
	}

	return result;
}


pixCompare computeSSD(pixCompare feature, R2Image &img, R2Image &other) {

	double soFarSSD = -1.0;
	double tempSSD = 0.0;

	int x_coord = feature.x;
	int y_coord = feature.y;

	int searchWindowWidthSpan = (other.Width() / 5) / 2;
	int searchWindowHeightSpan = (other.Height() / 5) / 2;

	pixCompare result;
	int match_x;
	int match_y;

	int featureSpan = 10;

	int dx, dy;

	for (dx = -1 * searchWindowWidthSpan; dx <= searchWindowWidthSpan; dx++) {
		for (dy = -1 * searchWindowHeightSpan; dy <= searchWindowHeightSpan; dy++) {

			if (x_coord + dx - featureSpan >= 0 && x_coord + dx + featureSpan < other.Width() &&
				y_coord + dy - featureSpan >= 0 && y_coord + dy + featureSpan < other.Height()) {

				tempSSD = singleSSD(x_coord, y_coord, x_coord + dx, y_coord + dy, img, other);

				if (tempSSD < soFarSSD || soFarSSD < 0.0) {
					match_x = dx + x_coord;
					match_y = dy + y_coord;
					soFarSSD = tempSSD;
				}

				tempSSD = 0.0;
			}
		}
	}

	// save the details of the best feature matching candidate and return it
	result.x = match_x;
	result.y = match_y;
	result.pix = other.Pixel(match_x, match_y);

	return result;
}



void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// computes Harris
	R2Image copy(*this);

	R2Pixel offset(0.5, 0.5, 0.5, 1);

	R2Image sobel_X(*this);
	sobel_X.SobelX();
	R2Image sobel_Y(*this);
	sobel_Y.SobelY();
	R2Image sobel_XY(*this);


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			sobel_XY.Pixel(i, j) = sobel_X.Pixel(i, j) * sobel_Y.Pixel(i, j);
			sobel_X.Pixel(i, j) = sobel_X.Pixel(i, j) * sobel_X.Pixel(i, j);
			sobel_Y.Pixel(i, j) = sobel_Y.Pixel(i, j) * sobel_Y.Pixel(i, j);
		}
	}

	sobel_X.Blur(2.0);
	sobel_Y.Blur(2.0);
	sobel_XY.Blur(2.0);

	R2Pixel tempPix;
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			tempPix = sobel_X.Pixel(x, y) * sobel_Y.Pixel(x, y) - sobel_XY.Pixel(x, y) * sobel_XY.Pixel(x, y) -
				0.04 * ((sobel_X.Pixel(x, y) + sobel_Y.Pixel(x, y)) * (sobel_X.Pixel(x, y) + sobel_Y.Pixel(x, y))) + offset;

			copy.Pixel(x, y) = tempPix;
		}
	}

	priority_queue<pixCompare, vector<pixCompare>, Compare> pq;
	pq = MakePriorityQueue(copy);

	// get the coordinates of the top features
	pixCompare* features = FeatureDetectionNew(copy, *this, pq);

	pixCompare otherFeat;

	pixCompare* matchedFeatures = new pixCompare[150];


	for (int i = 0; i < 150; i++) {
		printf("dealing with feature # %d\n", i);
		otherFeat = computeSSD(features[i], *this, *otherImage);
		matchedFeatures[i] = otherFeat;
		printf("previous feature coord: x = %d :: y = %d\n", features[i].x, features[i].y);
		printf("matched feature coords: x = %d :: y = %d\n\n", otherFeat.x, otherFeat.y);
		line(features[i].x, otherFeat.x, features[i].y, otherFeat.y, 1.0, 0.0, 0.0);
	}

	// rejecting bad matches
	int randVectorIndex;
	pair <int, int> randVector;
	pair <int, int> otherVector;
	int treshold = 5;
	int maxInliers = 0;
	int currentInliers = 0;


	for (int j = 0; j < 10; j++) {
		currentInliers = 0;
		randVectorIndex = rand() % 150;
		randVector.first = abs(matchedFeatures[randVectorIndex].x - features[randVectorIndex].x);
		randVector.second = abs(matchedFeatures[randVectorIndex].y - features[randVectorIndex].y);


		for (int k = 0; k < 150; k++) {
			otherVector.first = abs(matchedFeatures[k].x - features[k].x);
			otherVector.second = abs(matchedFeatures[k].y - features[k].y);
			if (sqrt((randVector.first - otherVector.first) * (randVector.first - otherVector.first) +
				(randVector.second - otherVector.second) * (randVector.second - otherVector.second)) < treshold) {
				currentInliers++;
				line(matchedFeatures[k].x, features[k].x, matchedFeatures[k].y, features[k].y, 0.0, 1.0, 0.0);
			}
		}

		if (currentInliers > maxInliers) {
			maxInliers = currentInliers;
			for (int k = 0; k < 150; k++) {
				otherVector.first = abs(matchedFeatures[k].x - features[k].x);
				otherVector.second = abs(matchedFeatures[k].y - features[k].y);
				if (sqrt((randVector.first - otherVector.first) * (randVector.first - otherVector.first) +
					(randVector.second - otherVector.second) * (randVector.second - otherVector.second)) < treshold) {
					currentInliers++;
					line(matchedFeatures[k].x, features[k].x, matchedFeatures[k].y, features[k].y, 0.0, 1.0, 0.0);
				}
				else {
					line(matchedFeatures[k].x, features[k].x, matchedFeatures[k].y, features[k].y, 1.0, 0.0, 0.0);
				}
			}
		}

	}

}



// Given a 'from' and 'to' point of a transformation, compute the two equations
//   and add them to the matrix A
void computePairEquations(pair <int, int> from, pair <int, int> to, double **A, int index) {

	int Ax = from.first;
	int Ay = from.second;
	int Bx = to.first;
	int By = to.second;

	double firstEq[9] = { Ax, Ay, 1, 0, 0, 0, (-1) * (Bx * Ax), (-1) * (Bx * Ay), (-1) * Bx };
	double secondEq[9] = { 0, 0, 0, Ax, Ay, 1, (-1) * (By * Ax), (-1) * (By * Ay), (-1) * By };


	for (int i = 0; i < 9; i++) {
		for (int j = 1; j < 10; j++) {
			A[index][j] = firstEq[j - 1];
			A[index+1][j] = secondEq[j - 1];
		}
	}

}

// Using the computePairEquations function, create the full matrix A with all of the equations
void computeEquationMatrix(pair<int, int> *from, pair<int, int> *to, double **equations) {

	int index = 1;
	for (int i = 0; i < 4; i++) {
		pair<int, int> from_point = from[i];
		pair<int, int> to_point = to[i];

		computePairEquations(from_point, to_point, equations, index);
		index = index + 2;
	}

}


// Compute the SVD of matrix A, find the singular values, and extract the nullspace
double** computeHomography(pair<int, int> *from, pair<int, int> *to) {
	
	double **equations;
	equations = new double *[9];
	for (int i = 0; i < 9; i++) {
		equations[i] = new double[10];
		for (int j = 0; j < 10; j++) {
			equations[i][j] = 0.0;
		}
	}

	computeEquationMatrix(from, to, equations);

	double singularValues[10]; // 1..9
	double** nullspaceMatrix = dmatrix(1, 9, 1, 9);


	svdcmp(equations, 8, 9, singularValues, nullspaceMatrix);

	// find the smallest singular value:
	int smallestIndex = 1;
	for (int i = 2; i< 10; i++) if (singularValues[i]<singularValues[smallestIndex]) smallestIndex = i;

	double **homographyMatrix;
	homographyMatrix = new double *[3];
	for (int i = 0; i < 3; i++) {
		homographyMatrix[i] = new double[3];
		for (int j = 0; j < 3; j++) {
			homographyMatrix[i][j] = 0.0;
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			homographyMatrix[i][j] = nullspaceMatrix[(3 * i) + j + 1][smallestIndex];
		}
	}

	return homographyMatrix;
}


void R2Image::
testHomography() {
	pair<int, int> from[4] = { make_pair(0, 0), make_pair(1, 0), make_pair(1, 1), make_pair(0, 1) };
	pair<int, int> to[4] = { make_pair(1, 0), make_pair(2, 0), make_pair(2, 1), make_pair(1, 1) };
	double** homographyMatrix = computeHomography(from, to);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f ", homographyMatrix[i][j]);
		}
		printf("\n");
	}
}


// RANSAC Algorithme using DLT 
boolean * R2Image::
dltRansac(R2Image * otherImage) {

	// apply Harris
	R2Image copy(*this);

	R2Pixel offset(0.5, 0.5, 0.5, 1);

	R2Image sobel_X(*this);
	sobel_X.SobelX();
	R2Image sobel_Y(*this);
	sobel_Y.SobelY();
	R2Image sobel_XY(*this);


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			sobel_XY.Pixel(i, j) = sobel_X.Pixel(i, j) * sobel_Y.Pixel(i, j);
			sobel_X.Pixel(i, j) = sobel_X.Pixel(i, j) * sobel_X.Pixel(i, j);
			sobel_Y.Pixel(i, j) = sobel_Y.Pixel(i, j) * sobel_Y.Pixel(i, j);
		}
	}

	sobel_X.Blur(2.0);
	sobel_Y.Blur(2.0);
	sobel_XY.Blur(2.0);

	R2Pixel tempPix;
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			tempPix = sobel_X.Pixel(x, y) * sobel_Y.Pixel(x, y) - sobel_XY.Pixel(x, y) * sobel_XY.Pixel(x, y) -
				0.04 * ((sobel_X.Pixel(x, y) + sobel_Y.Pixel(x, y)) * (sobel_X.Pixel(x, y) + sobel_Y.Pixel(x, y))) + offset;

			copy.Pixel(x, y) = tempPix;
		}
	}

	priority_queue<pixCompare, vector<pixCompare>, Compare> pq;
	pq = MakePriorityQueue(copy);

	// get the coordinates of the top features
	pixCompare* features = FeatureDetectionNew(copy, *this, pq);

	pixCompare* matchedFeatures = new pixCompare[150];

	// track the features
	for (int i = 0; i < 150; i++) {
		matchedFeatures[i] = computeSSD(features[i], *this, *otherImage);
	}

	double** homographyMatrix;
	static double bestHomographyMatrix[9];
	static boolean mostInliers[150];
	boolean tempInliers[150];

	double x, y, w;
	int maxInliers = 0;
	int threshold = 5;
	int currentInliers = 0;
	pair <double, double> initialVector;
	pair <double, double> homographyVector;

	int randPoint;
	pair<int, int> from_points[4];
	pair<int, int> to_points[4];

	double distance;

	// iterrate 1000 times
	for (int a = 0; a < 1000; a++) {
		currentInliers = 0;
		
		// compute the homography matrix H
		for (int p = 0; p < 4; p++) {
			randPoint = rand() % 150;
			from_points[p] = make_pair(features[randPoint].x, features[randPoint].y);
			to_points[p] = make_pair(matchedFeatures[randPoint].x, matchedFeatures[randPoint].y);
		}
		
		homographyMatrix = computeHomography(from_points, to_points);

		// using the computed matrix H, check other correspndences
		for (int i = 0; i < 150; i++) {
			distance = 0.0;
			x = homographyMatrix[0][0] * features[i].x +
				homographyMatrix[0][1] * features[i].y +
				homographyMatrix[0][2] * 1;
			y = homographyMatrix[1][0] * features[i].x +
				homographyMatrix[1][1] * features[i].y +
				homographyMatrix[1][2] * 1;
			w = homographyMatrix[2][0] * features[i].x +
				homographyMatrix[2][1] * features[i].y +
				homographyMatrix[2][2] * 1;
			x = x / w;
			y = y / w;

			initialVector.first = abs( matchedFeatures[i].x - features[i].x);
			initialVector.second = abs( matchedFeatures[i].y - features[i].y);

			homographyVector.first = abs(x - features[i].x);
			homographyVector.second = abs(y - features[i].y);

			distance = sqrt(((initialVector.first - homographyVector.first) * (initialVector.first - homographyVector.first)) +
				((initialVector.second - homographyVector.second) * (initialVector.second - homographyVector.second)));

			if (distance  < threshold) {
				currentInliers++;
				tempInliers[i] = true;
			}
			else {
				tempInliers[i] = false;
			}
		}

		if (currentInliers > maxInliers) {
			maxInliers = currentInliers;

			for (int m = 0; m < 150; m++) {
				mostInliers[m] = tempInliers[m];
			}
		}
	}

	for (int i = 0; i < 150; i++) {
		if (mostInliers[i]) {
			line(matchedFeatures[i].x, features[i].x, matchedFeatures[i].y, features[i].y, 0.0, 1.0, 0.0);
		}
		else {
			line(matchedFeatures[i].x, features[i].x, matchedFeatures[i].y, features[i].y, 1.0, 0.0, 0.0);
		}
	}

	return mostInliers;
}



void R2Image::
initialBlendOtherImageHomography(R2Image * otherImage)
{
	
	R2Image warped(width, height);
	double x_x, y_y, u;
	int warped_x, warped_y;

	boolean * mostInliers;
	R2Image copy(*this);
	
	mostInliers = dltRansac(otherImage);

	/*
	for (int w = 0; w < width; w++) {
		for (int h = 0; h < height; h++) {
			x_x = bestHomographyMatrix[0] * w +	bestHomographyMatrix[1] * h + bestHomographyMatrix[2] * 1;
			y_y = bestHomographyMatrix[3] * w +	bestHomographyMatrix[4] * h + bestHomographyMatrix[5] * 1;
			u = bestHomographyMatrix[6] * w + bestHomographyMatrix[7] * h + bestHomographyMatrix[8] * 1;
			
			warped_x = x_x / u;
			warped_x = min(warped_x, width-1);
			warped_x = max(0, warped_x);

			warped_y = y_y / u;
			warped_y = min(warped_y, height-1);
			warped_y = max(0, warped_y);

			warped.Pixel(w, h) = otherImage->Pixel(warped_x, warped_y);
		}
	}

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) = (warped.Pixel(i, j) + copy.Pixel(i, j))/2;
		}
	}*/
}


// Improved pixel resampling using bilinear pixel interpolation
void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	/*
	R2Image warped(width, height);
	R2Image copy(*this);

	R2Pixel p1, p2, p3, p4;
	R2Pixel warpedPix;

	double x, y, u;
	double dx, dy, dx1, dy1;

	int warped_x, warped_y;
	int w1, w2, w3, w4;
	int warp_r, warp_g, warp_b, warp_a;

	double * bestHomographyMatrix;

	// compute the best homography matrix
	bestHomographyMatrix = dltRansac(otherImage);


	for (int w = 0; w < width; w++) {
		for (int h = 0; h < height; h++) {
			x = bestHomographyMatrix[0] * w + bestHomographyMatrix[1] * h + bestHomographyMatrix[2] * 1;
			y = bestHomographyMatrix[3] * w + bestHomographyMatrix[4] * h + bestHomographyMatrix[5] * 1;
			u = bestHomographyMatrix[6] * w + bestHomographyMatrix[7] * h + bestHomographyMatrix[8] * 1;

			x = x / u;
			y = y / u;

			// Determining the position of the pixel
			warped_x = (int)x;
			warped_y = (int)y;

			// Getting the four neighboring pixels
			p1 = Pixel(warped_x, warped_y);
			p2 = Pixel(warped_x + 1, warped_y);
			p3 = Pixel(warped_x, warped_y + 1);
			p4 = Pixel(warped_x + 1, warped_y + 1);

			// Calculating the weights
			dx = x - warped_x;
			dy = y - warped_y;
			dx1 = 1.0 - dx;
			dy1 = 1.0 - dy;

			w1 = dx1 * dy1 * 256.0;
			w2 = dx * dy1 * 256.0;
			w3 = dx1 * dy * 256.0;
			w4 = dx * dy * 256.0;

			// Calculating the weighted sum of pixels
			warp_r = p1.Red * w1 + p2.Red * w2 + p3.Red * w3 + p4.Red * w4;
			warp_g = p1.Green * w1 + p2.Green * w2 + p3.Green * w3 + p4.Green * w4;
			warp_b = p1.Blue * w1 + p2.Blue * w2 + p3.Blue * w3 + p4.Blue * w4;
			warp_a = p1.Alpha * w1 + p2.Alpha * w2 + p3.Alpha * w3 + p4.Alpha * w4;

			// Assigning the value to the new pixel
			warpedPix.Red = warp_r;
			warpedPix.Green = warp_g;
			warpedPix.Blue = warp_b;
			warpedPix.Alpha = warp_a;

			warped.Pixel(w, h) = warpedPix;
		}
	}

	// blend images
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) = (warped.Pixel(i, j) + copy.Pixel(i, j)) / 2;
		}
	}*/
}


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
	// Initialize everything
	if (pixels) { delete[] pixels; pixels = NULL; }
	npixels = width = height = 0;

	// Parse input filename extension
	char *input_extension;
	if (!(input_extension = (char*)strrchr(filename, '.'))) {
		fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
		return 0;
	}

	// Read file of appropriate type
	if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
	else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
	else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
	else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

	// Should never get here
	fprintf(stderr, "Unrecognized image file extension");
	return 0;
}



int R2Image::
Write(const char *filename) const
{
	// Parse input filename extension
	char *input_extension;
	if (!(input_extension = (char*)strrchr(filename, '.'))) {
		fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
		return 0;
	}

	// Write file of appropriate type
	if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
	else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
	else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
	else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

	// Should never get here
	fprintf(stderr, "Unrecognized image file extension");
	return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
	unsigned short int bfType;
	unsigned int bfSize;
	unsigned short int bfReserved1;
	unsigned short int bfReserved2;
	unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
	unsigned int biSize;
	int biWidth;
	int biHeight;
	unsigned short int biPlanes;
	unsigned short int biBitCount;
	unsigned int biCompression;
	unsigned int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	unsigned int biClrUsed;
	unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
	unsigned char rgbtBlue;
	unsigned char rgbtGreen;
	unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
	unsigned char rgbBlue;
	unsigned char rgbGreen;
	unsigned char rgbRed;
	unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
	// Read a unsigned short int from a file in little endian format 
	unsigned short int lsb, msb;
	lsb = getc(fp);
	msb = getc(fp);
	return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
	// Write a unsigned short int to a file in little endian format
	unsigned char lsb = (unsigned char)(x & 0x00FF); putc(lsb, fp);
	unsigned char msb = (unsigned char)(x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
	// Read a unsigned int word from a file in little endian format 
	unsigned int b1 = getc(fp);
	unsigned int b2 = getc(fp);
	unsigned int b3 = getc(fp);
	unsigned int b4 = getc(fp);
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
	// Write a unsigned int to a file in little endian format 
	unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
	unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
	unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
	unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
	// Read a int word from a file in little endian format 
	int b1 = getc(fp);
	int b2 = getc(fp);
	int b3 = getc(fp);
	int b4 = getc(fp);
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
	// Write a int to a file in little endian format 
	char b1 = (x & 0x000000FF); putc(b1, fp);
	char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
	char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
	char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	/* Read file header */
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = WordReadLE(fp);
	bmfh.bfSize = DWordReadLE(fp);
	bmfh.bfReserved1 = WordReadLE(fp);
	bmfh.bfReserved2 = WordReadLE(fp);
	bmfh.bfOffBits = DWordReadLE(fp);

	/* Check file header */
	assert(bmfh.bfType == BMP_BF_TYPE);
	/* ignore bmfh.bfSize */
	/* ignore bmfh.bfReserved1 */
	/* ignore bmfh.bfReserved2 */
	assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

	/* Read info header */
	BITMAPINFOHEADER bmih;
	bmih.biSize = DWordReadLE(fp);
	bmih.biWidth = LongReadLE(fp);
	bmih.biHeight = LongReadLE(fp);
	bmih.biPlanes = WordReadLE(fp);
	bmih.biBitCount = WordReadLE(fp);
	bmih.biCompression = DWordReadLE(fp);
	bmih.biSizeImage = DWordReadLE(fp);
	bmih.biXPelsPerMeter = LongReadLE(fp);
	bmih.biYPelsPerMeter = LongReadLE(fp);
	bmih.biClrUsed = DWordReadLE(fp);
	bmih.biClrImportant = DWordReadLE(fp);

	// Check info header 
	assert(bmih.biSize == BMP_BI_SIZE);
	assert(bmih.biWidth > 0);
	assert(bmih.biHeight > 0);
	assert(bmih.biPlanes == 1);
	assert(bmih.biBitCount == 24);  /* RGB */
	assert(bmih.biCompression == BI_RGB);   /* RGB */
	int lineLength = bmih.biWidth * 3;  /* RGB */
	if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
	assert(bmih.biSizeImage == (unsigned int)lineLength * (unsigned int)bmih.biHeight);

	// Assign width, height, and number of pixels
	width = bmih.biWidth;
	height = bmih.biHeight;
	npixels = width * height;

	// Allocate unsigned char buffer for reading pixels
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = bmih.biSizeImage;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Read buffer 
	fseek(fp, (long)bmfh.bfOffBits, SEEK_SET);
	if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
		fprintf(stderr, "Error while reading BMP file %s", filename);
		return 0;
	}

	// Close file
	fclose(fp);

	// Allocate pixels for image
	pixels = new R2Pixel[width * height];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double b = (double)*(p++) / 255;
			double g = (double)*(p++) / 255;
			double r = (double)*(p++) / 255;
			R2Pixel pixel(r, g, b, 1);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return success
	return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
	// Open file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Compute number of bytes in row
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

	// Write file header 
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = BMP_BF_TYPE;
	bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfOffBits = BMP_BF_OFF_BITS;
	WordWriteLE(bmfh.bfType, fp);
	DWordWriteLE(bmfh.bfSize, fp);
	WordWriteLE(bmfh.bfReserved1, fp);
	WordWriteLE(bmfh.bfReserved2, fp);
	DWordWriteLE(bmfh.bfOffBits, fp);

	// Write info header 
	BITMAPINFOHEADER bmih;
	bmih.biSize = BMP_BI_SIZE;
	bmih.biWidth = width;
	bmih.biHeight = height;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;       /* RGB */
	bmih.biCompression = BI_RGB;    /* RGB */
	bmih.biSizeImage = rowsize * (unsigned int)bmih.biHeight;  /* RGB */
	bmih.biXPelsPerMeter = 2925;
	bmih.biYPelsPerMeter = 2925;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;
	DWordWriteLE(bmih.biSize, fp);
	LongWriteLE(bmih.biWidth, fp);
	LongWriteLE(bmih.biHeight, fp);
	WordWriteLE(bmih.biPlanes, fp);
	WordWriteLE(bmih.biBitCount, fp);
	DWordWriteLE(bmih.biCompression, fp);
	DWordWriteLE(bmih.biSizeImage, fp);
	LongWriteLE(bmih.biXPelsPerMeter, fp);
	LongWriteLE(bmih.biYPelsPerMeter, fp);
	DWordWriteLE(bmih.biClrUsed, fp);
	DWordWriteLE(bmih.biClrImportant, fp);

	// Write image, swapping blue and red in each pixel
	int pad = rowsize - width * 3;
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			const R2Pixel& pixel = (*this)[i][j];
			double r = 255.0 * pixel.Red();
			double g = 255.0 * pixel.Green();
			double b = 255.0 * pixel.Blue();
			if (r >= 255) r = 255;
			if (g >= 255) g = 255;
			if (b >= 255) b = 255;
			fputc((unsigned char)b, fp);
			fputc((unsigned char)g, fp);
			fputc((unsigned char)r, fp);
		}

		// Pad row
		for (int i = 0; i < pad; i++) fputc(0, fp);
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Read PPM file magic identifier
	char buffer[128];
	if (!fgets(buffer, 128, fp)) {
		fprintf(stderr, "Unable to read magic id in PPM file");
		fclose(fp);
		return 0;
	}

	// skip comments
	int c = getc(fp);
	while (c == '#') {
		while (c != '\n') c = getc(fp);
		c = getc(fp);
	}
	ungetc(c, fp);

	// Read width and height
	if (fscanf(fp, "%d%d", &width, &height) != 2) {
		fprintf(stderr, "Unable to read width and height in PPM file");
		fclose(fp);
		return 0;
	}

	// Read max value
	double max_value;
	if (fscanf(fp, "%lf", &max_value) != 1) {
		fprintf(stderr, "Unable to read max_value in PPM file");
		fclose(fp);
		return 0;
	}

	// Allocate image pixels
	pixels = new R2Pixel[width * height];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for PPM file");
		fclose(fp);
		return 0;
	}

	// Check if raw or ascii file
	if (!strcmp(buffer, "P6\n")) {
		// Read up to one character of whitespace (\n) after max_value
		int c = getc(fp);
		if (!isspace(c)) putc(c, fp);

		// Read raw image data 
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				double r = (double)getc(fp) / max_value;
				double g = (double)getc(fp) / max_value;
				double b = (double)getc(fp) / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}
	else {
		// Read asci image data 
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				// Read pixel values
				int red, green, blue;
				if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
					fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
					fclose(fp);
					return 0;
				}

				// Assign pixel values
				double r = (double)red / max_value;
				double g = (double)green / max_value;
				double b = (double)blue / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
	// Check type
	if (ascii) {
		// Open file
		FILE *fp = fopen(filename, "w");
		if (!fp) {
			fprintf(stderr, "Unable to open image file: %s", filename);
			return 0;
		}

		// Print PPM image file 
		// First ppm pixel is top-left, so write in opposite scan-line order
		fprintf(fp, "P3\n");
		fprintf(fp, "%d %d\n", width, height);
		fprintf(fp, "255\n");
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int)(255 * p.Red());
				int g = (int)(255 * p.Green());
				int b = (int)(255 * p.Blue());
				fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
				if (((i + 1) % 4) == 0) fprintf(fp, "\n");
			}
			if ((width % 4) != 0) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");

		// Close file
		fclose(fp);
	}
	else {
		// Open file
		FILE *fp = fopen(filename, "wb");
		if (!fp) {
			fprintf(stderr, "Unable to open image file: %s", filename);
			return 0;
		}

		// Print PPM image file 
		// First ppm pixel is top-left, so write in opposite scan-line order
		fprintf(fp, "P6\n");
		fprintf(fp, "%d %d\n", width, height);
		fprintf(fp, "255\n");
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int)(255 * p.Red());
				int g = (int)(255 * p.Green());
				int b = (int)(255 * p.Blue());
				fprintf(fp, "%c%c%c", r, g, b);
			}
		}

		// Close file
		fclose(fp);
	}

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
};
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Initialize decompression info
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, fp);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	// Remember image attributes
	width = cinfo.output_width;
	height = cinfo.output_height;
	npixels = width * height;
	int ncomponents = cinfo.output_components;

	// Allocate pixels for image
	pixels = new R2Pixel[npixels];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Allocate unsigned char buffer for reading image
	int rowsize = ncomponents * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
		fclose(fp);
		return 0;
	}

	// Read scan lines 
	// First jpeg pixel is top-left, so read pixels in opposite scan-line order
	while (cinfo.output_scanline < cinfo.output_height) {
		int scanline = cinfo.output_height - cinfo.output_scanline - 1;
		unsigned char *row_pointer = &buffer[scanline * rowsize];
		jpeg_read_scanlines(&cinfo, &row_pointer, 1);
	}

	// Free everything
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	// Close file
	fclose(fp);

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double r, g, b, a;
			if (ncomponents == 1) {
				r = g = b = (double)*(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 1) {
				r = g = b = (double)*(p++) / 255;
				a = 1;
				p++;
			}
			else if (ncomponents == 3) {
				r = (double)*(p++) / 255;
				g = (double)*(p++) / 255;
				b = (double)*(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 4) {
				r = (double)*(p++) / 255;
				g = (double)*(p++) / 255;
				b = (double)*(p++) / 255;
				a = (double)*(p++) / 255;
			}
			else {
				fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
				return 0;
			}
			R2Pixel pixel(r, g, b, a);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return success
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
	// Open file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Initialize compression info
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width = width; 	/* image width and height, in pixels */
	cinfo.image_height = height;
	cinfo.input_components = 3;		/* # of color components per pixel */
	cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
	cinfo.dct_method = JDCT_ISLOW;
	jpeg_set_defaults(&cinfo);
	cinfo.optimize_coding = TRUE;
	jpeg_set_quality(&cinfo, 95, TRUE);
	jpeg_start_compress(&cinfo, TRUE);

	// Allocate unsigned char buffer for reading image
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
		fclose(fp);
		return 0;
	}

	// Fill buffer with pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			const R2Pixel& pixel = (*this)[i][j];
			int r = (int)(255 * pixel.Red());
			int g = (int)(255 * pixel.Green());
			int b = (int)(255 * pixel.Blue());
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			*(p++) = r;
			*(p++) = g;
			*(p++) = b;
		}
	}



	// Output scan lines
	// First jpeg pixel is top-left, so write in opposite scan-line order
	while (cinfo.next_scanline < cinfo.image_height) {
		int scanline = cinfo.image_height - cinfo.next_scanline - 1;
		unsigned char *row_pointer = &buffer[scanline * rowsize];
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}

	// Free everything
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	// Close file
	fclose(fp);

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return number of bytes written
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}






