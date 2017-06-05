
// Computer Vision for Digital Post-Production
// Lecturer: Gergely Vass - vassg@vassg.hu
//
// Skeleton Code for programming assigments
//
// Code originally from Thomas Funkhouser
// main.c
// original by Wagner Correa, 1999
// modified by Robert Osada, 2000
// modified by Renato Werneck, 2003
// modified by Jason Lawrence, 2004
// modified by Jason Lawrence, 2005
// modified by Forrester Cole, 2006
// modified by Tom Funkhouser, 2007
// modified by Chris DeCoro, 2007
//



// Include files

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <sstream>
#include <time.h>
#include <map>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"

using namespace std;


// Program arguments

static char options[] =
"  -help\n"
"  -svdTest\n"
"  -sobelX\n"
"  -sobelY\n"
"  -log\n"
"  -harris <real:sigma>\n"
"  -saturation <real:factor>\n"
"  -brightness <real:factor>\n"
"  -blur <real:sigma>\n"
"  -sharpen \n"
"  -rejectBadMatches \n"
"  -dltRansac \n"
"  -computeHomography \n"
"  -matchTranslation <file:other_image>\n"
"  -matchHomography <file:other_image>\n";


static void
ShowUsage(void)
{
	// Print usage message and exit
	fprintf(stderr, "Usage: imgpro input_image output_image [  -option [arg ...] ...]\n");
	fprintf(stderr, options);
	exit(EXIT_FAILURE);
}



static void
CheckOption(char *option, int argc, int minargc)
{
	// Check if there are enough remaining arguments for option
	if (argc < minargc) {
		fprintf(stderr, "Too few arguments for %s\n", option);
		ShowUsage();
		exit(-1);
	}
}



static int
ReadCorrespondences(char *filename, R2Segment *&source_segments, R2Segment *&target_segments, int& nsegments)
{
	// Open file
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		fprintf(stderr, "Unable to open correspondences file %s\n", filename);
		exit(-1);
	}

	// Read number of segments
	if (fscanf(fp, "%d", &nsegments) != 1) {
		fprintf(stderr, "Unable to read correspondences file %s\n", filename);
		exit(-1);
	}

	// Allocate arrays for segments
	source_segments = new R2Segment[nsegments];
	target_segments = new R2Segment[nsegments];
	if (!source_segments || !target_segments) {
		fprintf(stderr, "Unable to allocate correspondence segments for %s\n", filename);
		exit(-1);
	}

	// Read segments
	for (int i = 0; i < nsegments; i++) {

		// Read source segment
		double sx1, sy1, sx2, sy2;
		if (fscanf(fp, "%lf%lf%lf%lf", &sx1, &sy1, &sx2, &sy2) != 4) {
			fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
			exit(-1);
		}

		// Read target segment
		double tx1, ty1, tx2, ty2;
		if (fscanf(fp, "%lf%lf%lf%lf", &tx1, &ty1, &tx2, &ty2) != 4) {
			fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
			exit(-1);
		}

		// Add segments to list
		source_segments[i] = R2Segment(sx1, sy1, sx2, sy2);
		target_segments[i] = R2Segment(tx1, ty1, tx2, ty2);
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



int
main(int argc, char **argv) {

	srand(time(NULL));
	// Look for help
	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-help")) {
			ShowUsage();
		}
		if (!strcmp(argv[i], "-svdTest")) {
			R2Image *image = new R2Image();
			image->svdTest();
			return 0;
		}
	}

	// Read input and output image filenames
	if (argc < 3)  ShowUsage();
	argv++, argc--; // First argument is program name
	char *input_image_name = *argv; argv++, argc--;
	char *output_image_name = *argv; argv++, argc--;

	// Allocate image
	R2Image *image = new R2Image();
	if (!image) {
		fprintf(stderr, "Unable to allocate image\n");
		exit(-1);
	}

	// Read input image
	/*if (!image->Read(input_image_name)) {
	fprintf(stderr, "Unable to read image from %s\n", input_image_name);
	exit(-1);
	}*/

	// Initialize sampling method
	int sampling_method = R2_IMAGE_POINT_SAMPLING;

	// Parse arguments and perform operations
	while (argc > 0) {
		if (!strcmp(*argv, "-brightness")) {
			CheckOption(*argv, argc, 2);
			double factor = atof(argv[1]);
			argv += 2, argc -= 2;
			image->Brighten(factor);
		}
		else if (!strcmp(*argv, "-sobelX")) {
			argv++, argc--;
			image->SobelX();
		}
		else if (!strcmp(*argv, "-sobelY")) {
			argv++, argc--;
			image->SobelY();
		}
		else if (!strcmp(*argv, "-log")) {
			argv++, argc--;
			image->LoG();
		}
		else if (!strcmp(*argv, "-saturation")) {
			CheckOption(*argv, argc, 2);
			double factor = atof(argv[1]);
			argv += 2, argc -= 2;
			image->ChangeSaturation(factor);
		}
		else if (!strcmp(*argv, "-harris")) {
			CheckOption(*argv, argc, 2);
			double sigma = atof(argv[1]);
			argv += 2, argc -= 2;
			image->Harris(sigma);
		}
		else if (!strcmp(*argv, "-blur")) {
			CheckOption(*argv, argc, 2);
			double sigma = atof(argv[1]);
			argv += 2, argc -= 2;
			image->Blur(sigma);
		}
		else if (!strcmp(*argv, "-sharpen")) {
			argv++, argc--;
			image->Sharpen();
		}
		else if (!strcmp(*argv, "-testHomography")) {
			argv++, argc--;
			image->testHomography();
		}
		else if (!strcmp(*argv, "-dltRansac")) {
			CheckOption(*argv, argc, 2);
			R2Image *other_image = new R2Image(argv[1]);
			argv += 2, argc -= 2;
			image->dltRansac(other_image);
			delete other_image;
		}
		else if (!strcmp(*argv, "-matchTranslation")) {
			CheckOption(*argv, argc, 2);
			R2Image *other_image = new R2Image(argv[1]);
			argv += 2, argc -= 2;
			image->blendOtherImageTranslated(other_image);
			delete other_image;
		}
		else if (!strcmp(*argv, "-matchHomography")) {
			CheckOption(*argv, argc, 2);
			R2Image *other_image = new R2Image(argv[1]);
			argv += 2, argc -= 2;
			image->blendOtherImageHomography(other_image);
			delete other_image;
		}

		else if (!strcmp(*argv, "-Final")) {
			int maxFrame = 360;
			const int count = 4;
			//R2Image* toFit = new R2Image(argv[1]);
			R2Image* image = new R2Image((string(input_image_name) + "115.jpg").c_str());
			R2Image* firstSeq = new R2Image((string(input_image_name) + "t.jpg").c_str());
			R2Image* mask = new R2Image("../../Input/seq1/mask.jpg");

			R2Image *other_image = 0;
			R2Image *copyMask = 0;
			R2Image *copyFirstseq = 0;

			int featureWidth = 60;
			int featureHeight = 60;

			pair<int, int> features[count];
			pair<int, int> matched[count];
			pair<int, int> first[count];

			first[0] = make_pair(655, 830);
			first[1] = make_pair(577, 692);
			first[2] = make_pair(732, 699);
			first[3] = make_pair(665, 609);
			//first[4] = make_pair(661, 725);

			features[0] = first[0];
			features[1] = first[1];
			features[2] = first[2];
			features[3] = first[3];
			//features[4] = first[4];

			for (int frame = 115; frame <= maxFrame; frame++) {
				cout << "\t" << frame << "\n";
				string frameNum = to_string(frame);
				other_image = new R2Image((input_image_name + frameNum + ".jpg").c_str());

				for (int k = 0; k < count; k++) {
					matched[k] = image->finalTrack(other_image, features[k], featureWidth, featureHeight);
					features[k] = matched[k];
				}

				//double currE = 0, BigError = -1;
				//int bestToIgnore = 4;
				/*
				for (int ignore = 0; ignore < count; ignore++) {
				currE = R2Image::RankHomography(first, matched, ignore, count);
				if (BigError < 0 || currE > BigError) {
				BigError = currE;
				bestToIgnore = ignore;
				}
				}*/

				/*
				pair<int, int> pfeature[4], pmatched[4];
				int selected = 0;
				for (int s = 0; s<count; s++) {
				if (s != bestToIgnore) {
				pfeature[selected] = first[s];
				pmatched[selected] = matched[s];
				selected++;
				}
				}*/

				double** homographyM = R2Image::computeHomographyFinal(matched, first);

				//other_image->MarkFeatures(matched, featureHeight / 2, 5);
				//other_image->MarkFeatures(pmatched, featureHeight / 2, 4);

				copyMask = new R2Image(*mask);
				copyFirstseq = new R2Image(*firstSeq);

				//other_image->ApplyMatrix(homographyM);
				copyMask->ApplyMatrix(homographyM);
				copyFirstseq->ApplyMatrix(homographyM);

				ostringstream ostr2;
				ostr2 << output_image_name << frame - 114;


				copyFirstseq->Write((ostr2.str() + ".jpg").c_str());
				copyMask->Write((ostr2.str() + "mask.jpg").c_str());
				other_image->Write((ostr2.str() + "n.jpg").c_str());

				//other_image->BlendWithMask(copyMask, copyFirstseq);
				//other_image->Write((ostr2.str() + ".jpg").c_str());

				delete image;
				delete copyMask;
				delete copyFirstseq;
				image = other_image;

			}

			cout << "Done\n";

			//delete other_image;
			return EXIT_SUCCESS;
		}

		else if (!strcmp(*argv, "-blend")) {
			int maxFrame = 246;
			R2Image *original, *morphed, *mask;
			for (int frame = 1; frame <= maxFrame; frame++) {
				cout << frame << endl;
				original = new R2Image((string(input_image_name) + to_string(frame) + "n.jpg").c_str());
				morphed = new R2Image((string(input_image_name) + to_string(frame) + ".jpg").c_str());
				mask = new R2Image((string(input_image_name) + to_string(frame) + "mask.jpg").c_str());

				original->BlendWithMask(mask, morphed);

				original->Write((string(output_image_name) + to_string(frame) + "o.jpg").c_str());
				delete original;
				delete morphed;
				delete mask;
			}
			return EXIT_SUCCESS;
		}
		else {
			// Unrecognized program argument
			fprintf(stderr, "image: invalid option: %s\n", *argv);
			ShowUsage();
		}
	}

	// Write output image
	if (!image->Write(output_image_name)) {
		fprintf(stderr, "Unable to read image from %s\n", output_image_name);
		exit(-1);
	}

	// Delete image
	delete image;

	// Return success
	return EXIT_SUCCESS;
}
