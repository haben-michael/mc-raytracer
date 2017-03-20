// Source file for image class



// Include files 

#include "R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
//#include <math.h>
#include <cmath>



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
	pixels = new R2Pixel [ npixels ];
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
	pixels = new R2Pixel [ npixels ];
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
	pixels = new R2Pixel [ npixels ];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++) 
		pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
	// Free image pixels
	if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
	// Delete previous pixels
	if (pixels) { delete [] pixels; pixels = NULL; }

	// Reset width and height
	npixels = image.npixels;
	width = image.width;
	height = image.height;

	// Allocate new pixels
	pixels = new R2Pixel [ npixels ];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++) 
		pixels[i] = image.pixels[i];

	// Return image
	return *this;
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
		for (int j = 0;  j < height; j++) {
			Pixel(i,j) *= factor;
			Pixel(i,j).Clamp();
		}
	}
}

void R2Image::
AddNoise(double factor)
{
	// Add noise to an image.  The amount of noise is given by the factor
	// in the range [0.0..1.0].  0.0 adds no noise.  1.0 adds a lot of noise.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "AddNoise(%g) not implemented\n", factor);
	//srand((unsigned)time(0));
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++)
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++)
				if (double(rand())/((double)RAND_MAX+1) < factor) 
					Pixel(i,j).Components()[k]=(double)rand()/RAND_MAX;

}

void R2Image::
Speckle(double percentage)
{
	// Add speckled noise to an image.  The percentage of pixels that are speckled
	// is given by the parameter in the range [0.0..1.0].  0.0 changes no pixels,
	// 1.0 changes all the pixels to a random value (a different random value for each)
	// 0.5 will change half of the pixels.

	// MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
	fprintf(stderr, "Speckle(%g) not implemented\n", percentage);
}

void R2Image::
ChangeContrast(double factor)
{
	// Change the contrast of an image by interpolating between the image
	// and a constant gray image with the average luminance.
	// Interpolation reduces constrast, extrapolation boosts constrast,
	// and negative factors generate inverted images.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "ChangeContrast(%g) not implemented\n", factor);
	double aveLuminance = 0;
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++) 
			aveLuminance += Pixel(i,j).Luminance();
	aveLuminance /= this->NPixels();

	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++) 
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++) {
				Pixel(i,j).Components()[k] = (1-factor)*aveLuminance + factor*Pixel(i,j).Components()[k];
				Pixel(i,j).Clamp();
			}
}

void R2Image::
ChangeSaturation(double factor)
{
	// Changes the saturation of an image by interpolating between the
	// image and a gray level version of the image.  Interpolation
	// decreases saturation, extrapolation increases it, negative factors
	// preserve luminance  but invert the hue of the input image.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++) {
			double luminance=Pixel(i,j).Luminance();
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++) {
				//if (rand()/((double)RAND_MAX+1) < factor) 
				Pixel(i,j).Components()[k] = factor*Pixel(i,j).Components()[k] + (1-factor)*luminance;
				Pixel(i,j).Clamp();
			}
		}
}

void R2Image::
Threshold(double value)
{
	// Set all of the pixels with values less than value to 0.0

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Threshold(%g) not implemented\n", value);
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++) 
			if (Pixel(i,j).Luminance() < value)
				for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++) Pixel(i,j).Components()[k] = 0;

}

void R2Image::
ApplyGamma(double exponent)
{
	// Apply a gamma correction with exponent to each pixel

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Gamma(%g) not implemented\n", exponent);
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++)
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++) 
				Pixel(i,j).Components()[k] = pow(Pixel(i,j).Components()[k],exponent);
}

void R2Image::
BlackAndWhite(void)
{
	// Replace each pixel with its luminance value
	// Put this in each channel,  so the result is grayscale

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "BlackAndWhite() not implemented\n");
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++) {
			double luminance=Pixel(i,j).Luminance();
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++)
				Pixel(i,j).Components()[k] = luminance;
		}
}
//
//R2Image* toHSL(R2Image* rgbImage) {
//	typedef struct{double h,s,l;} ss[10];
//	assert(rgbImage!=NULL);
//	R2Image* hslImage = new R2Image(*rgbImage);
//	for (int i=0; i<hslImage->Width(); i++)
//		for (int j=0; j<hslImage->Height(); j++) {
//			double maxChan=max(hslImage->Pixel(i,j).Components()[0],
//				max(hslImage->Pixel(i,j).Components()[1],hslImage->Pixel(i,j).Components()[2]));
//			double minChan=min(hslImage->Pixel(i,j).Components()[0],
//				min(hslImage->Pixel(i,j).Components()[1],hslImage->Pixel(i,j).Components()[2]));
//			double l = .5*(maxChan-minChan);
//			double s = (l<=.5) ? (maxChan-minChan)/(2*l) : (maxChan-minChan)/(2-2*l);
//			//			if (maxChan==
//		}
//		return hslImage;
//}


double *toHSL(double r,double g, double b) {
	double *hsl = new double[3];
	double maxChan=max(r,max(g,b));
	double minChan=min(r,min(g,b));
	hsl[2] = .5*(maxChan+minChan);
	if (maxChan==minChan) hsl[1]=0;
	else
		hsl[1] = (hsl[2]<=.5) ? (maxChan-minChan)/(2*hsl[2]) : (maxChan-minChan)/(2-2*hsl[2]);
	if (maxChan==r)
		hsl[0] = fmod(60*(g-b)/(maxChan-minChan) + 360, 360);
	else
		if (maxChan==g)
			hsl[0] = 60*(b-r)/(maxChan-minChan) + 120;
		else
			if (maxChan==b)
				hsl[0] = 60*(r-g)/(maxChan-minChan) + 240;
			else hsl[0] = 0; //assume max=min
			return hsl;
}

double *toRGB(double h,double s,double l) {
	double *rgb = new double[3];
	double t_rgb[3];
	double q,p,h_k;

	q = (l<.5) ? l*(1+s) : l+s-l*s;
	p = 2*l-q;
	h_k = h/360;
	t_rgb[0] = h_k+(double)1/3;
	t_rgb[1] = h_k;
	t_rgb[2] = h_k-(double)1/3;
	double d=fmod(-.5,1);
	for (int k=0; k<3; k++) if (t_rgb[k]<0) t_rgb[k]+=1;
	for (int k=0; k<3; k++) if (t_rgb[k]>1) t_rgb[k]-=1;
	for (int k=0; k<3; k++) {
		if (t_rgb[k]<(double)1/6)
			rgb[k] = p+((q-p)*6*t_rgb[k]);
		else
			if (t_rgb[k]>=(double)1/6 && t_rgb[k]<(double)1/2)
				rgb[k] = q;
			else
				if (t_rgb[k]>=(double)1/2 && t_rgb[k]<(double)2/3)
					rgb[k] = p+((q-p)*6*((double)2/3-t_rgb[k]));
				else
					rgb[k] = p;
	}
	return rgb;
}

int* quantize(double* vals,int n_vals,int n_quantiles) {
	int *pdf = new int[n_quantiles];
	for (int i=0; i<n_quantiles; i++) pdf[i]=0;
	for (int i=0; i<n_vals; i++) {
		assert(vals[i]>=0 && vals[i]<=1);
		pdf[(int)floor((n_quantiles)*vals[i])]++;
	}
	return pdf;
}

void R2Image::
EqualizeLuminanceHistogram(void)
{
	// Increase the global contrast by converting to HLS color space, 
	// applying histogram equalization to the luminance channel, and 
	// converting back to RGB.
	// See http://en.wikipedia.org/wiki/Histogram_equalization

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "EqualizeLuminanceHistogram() not implemented\n");
	//	int n=50;
	//	double *a=new double[n];
	//	for (int i=0;i<n;i++) {a[i]=(double)rand()/RAND_MAX; printf("%.3f\t",a[i]);}
	//	int *x=quantize(a,n,10);
	//for (int i=0;i<10;i++) printf("%d\t",x[i]);

	//a=toHSL(1,1,1);
	//	a=toHSL(.5,1,.5);
	//	a=toHSL(0,0,.5);
	//a=toRGB(150.379,.050305,.92155);
	//a=toRGB(120,1,.75);
	//a=toRGB(240,1,.25);

	//create hsl array, with same structure as pixels array
	struct hsl_t {double h,s,l;};
	hsl_t *hslImage = new hsl_t[npixels];
	double *hsl;
	double r,g,b;

	//for (int i=0;i<10;i++) printf("%g,%g,%g\n",Pixel(1,i).Red(),Pixel(1,i).Green(),Pixel(1,i).Blue());

	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++) {
			r = Pixel(i,j).Components()[R2_IMAGE_RED_CHANNEL];
			g = Pixel(i,j).Components()[R2_IMAGE_GREEN_CHANNEL];
			b = Pixel(i,j).Components()[R2_IMAGE_BLUE_CHANNEL];
			hsl = toHSL(r,g,b);
			hslImage[i*height+j].h=hsl[0];
			hslImage[i*height+j].s=hsl[1];
			hslImage[i*height+j].l=hsl[2];
			double *rgb = toRGB(hsl[0],hsl[1],hsl[2]);
			if (!(abs(rgb[0]-r)<.001 && abs(rgb[1]-g)<.001 && abs(rgb[2]-b)<.001)) {
				printf("%g,%g,%g\n",r,g,b);
				printf("%g,%g,%g\n",rgb[0],rgb[1],rgb[2]);
				assert(false);

			}
		}

		//hsl conversion complete, now do equalization
		//first calculate pdf, quantizing luminance
		double h,s,l;
		int n_quantiles = 100;//min(npixels,80);
		double *luminances = new double[npixels];
		for (int i=0;i<npixels;i++) luminances[i] = hslImage[i].l;
		int *pdf = quantize(luminances,npixels,n_quantiles);

		//now calculate cdf
		int *cdf = new int[n_quantiles];
		cdf[0] = pdf[0];
		for (int i=1; i<n_quantiles; i++) cdf[i] = cdf[i-1]+pdf[i];
		//	printf("%d,%d",n_quantiles,npixels);
		/*for (int i=0;i<npixels;i++) {printf("%3.3f\t",hslImage[i].l);}
		printf("\n--------\n");
		for (int i=0;i<n_quantiles;i++) printf("%d\t",pdf[i]);
		printf("\n--------\n");
		for (int i=0;i<n_quantiles;i++) printf("%d\t",cdf[i]);
		printf("\n--------\n");*/

		//printf("TEST\n");
		//int lums[] = {52,55,61,66,70,61,64,73,
		//63,59,55,90,109,85,69,72,
		//62,59,68,113,144,104,66,73,
		//63,58,71,122,154,106,70,69,
		//67,61,68,104,126,88,68,70,
		//79,65,60,70,77,68,58,75,
		//85,71,64,59,55,61,65,83,
		//87,79,69,68,65,76,78,94};
		//int *pdf2 = new int[256];
		//for (int i=0;i<256;i++) pdf2[i]=0;
		//for (int i=0;i<64;i++) pdf2[lums[i]]++;
		//int *cdf2 = new int[256];
		//cdf2[0] = pdf2[0];
		//for (int i=1; i<256; i++) cdf2[i] = cdf2[i-1]+pdf2[i];
		//printf("\n--------\n");
		//for (int i=0;i<256;i++) printf("%d\t",cdf2[i]);
		//printf("\n--------\n");
		//int min_cdf2;
		//for (min_cdf2=0;cdf2[min_cdf2]==0 && min_cdf2<256;min_cdf2++);
		//for (int i=0; i<64; i++) {
		//	lums[i] = int((((double)cdf2[lums[i]]-cdf2[min_cdf2])/(64-cdf2[min_cdf2]))*(256-1)+.5);
		//}
		//for (int i=0;i<64;i++) printf("%d\t",lums[i]);

		//apply equalization formula
		int min_cdf;
		for (min_cdf=0;cdf[min_cdf]==0 && min_cdf<n_quantiles;min_cdf++);
		//printf("%d\n",min_cdf);
		for (int i=0,quantile; i<npixels; i++) {
			quantile = (int)floor((n_quantiles)*hslImage[i].l);
			//printf("%d:%.3f;",quantile,hslImage[i].l);
			hslImage[i].l = ((((double)cdf[quantile]-cdf[min_cdf])/(npixels-cdf[min_cdf]))*(n_quantiles-1))/n_quantiles;
			//printf("%.3f\n",hslImage[i].l);
		}
		//check on shape of new cdf
		for (int i=0;i<npixels;i++) luminances[i] = hslImage[i].l;
		int *pdf2 = quantize(luminances,npixels,n_quantiles);
		int *cdf2 = new int[n_quantiles];
		cdf2[0]=pdf2[0];
		for (int i=1; i<n_quantiles; i++) cdf2[i] = cdf2[i-1]+pdf2[i];
		for (int i=1; i<n_quantiles; i++) printf("%d\t",cdf2[i]);



		//now convert hsl back to rgb
		double *rgb;
		for (int i=0; i<width; i++)
			for (int j=0; j<height; j++) {
				h=hslImage[i*height+j].h;
				s=hslImage[i*height+j].s;
				l=hslImage[i*height+j].l;
				rgb = toRGB(h,s,l);
				for (int k=0; k<3; k++) Pixel(i,j).Components()[k]=rgb[k];
				//printf("-------\n");
				//			printf("%g,%g,%g\n",rgb[0],rgb[1],rgb[2]);

			}
			//for (int i=0;i<10;i++) printf("%g,%g,%g\n",Pixel(1,i).Red(),Pixel(1,i).Green(),Pixel(1,i).Blue());

			//for (int i=0; i<width; i++)
			//	for (int j=0; j<height; j++) {
			//		for (int k=0;k<3;k++) 
			//			Pixel(i,j).Components()[k]=0;
			//			Pixel(i,j).Components()[1]=100;

			//	}
			delete hslImage,luminances,pdf,cdf;

}


// Linear filtering ////////////////////////////////////////////////

class SumFilter;
class ProdFilter;

class Filter {
public: 
	double radius;
	virtual double f(double x) {return 0;}//don't want a complete virtual function so we can define +
	virtual double f(double x,double y) {return 0;}
	SumFilter operator+(Filter&);
	ProdFilter operator*(double);
};

class SumFilter : public Filter {
	Filter *filter1, *filter2;
public:
	SumFilter(Filter *filter1,Filter *filter2) {
		this->radius=max(filter1->radius,filter2->radius);
		this->filter1 = filter1;
		this->filter2 = filter2;
	}
	double f(double x) {
		//printf("%g",double(this->filter2->f(10)));
		return double(filter2->f(x)+filter1->f(x));
	}
};

class ProdFilter : public Filter {
	Filter *filter;
	double coef;
public:
	ProdFilter(Filter *filter,double coef) {
		this->radius=filter->radius;
		this->filter = filter;
		this->coef = coef;
	}
	double f(double x) {
		return double(coef*filter->f(x));
	}
};

SumFilter Filter::operator+(Filter &param)
{
	// add two filter functions
	return(SumFilter(this,&param));
}

ProdFilter Filter::operator*(double param)
{
	// add two filter functions
	return(ProdFilter(this,param));
}

class GaussianFilter : public Filter {
	double sigma;
	//int radius;
public:
	GaussianFilter(int radius,double sigma) {
		this->sigma = sigma;
		this->radius = radius;
	}
	GaussianFilter(double sigma) {
		this->sigma = sigma;
		this->radius = (int)floor(3.0*sigma+.5); //use 3*sigma rule
	}
	double f(double x) {
		return /*3.0*sigma/radius**/(double(1.0/sqrt(2*3.1416)/sigma)*exp(-x*x/sigma/sigma/2.0));
	}
	double f(double x,double y) {
		return this->f(x) * this->f(y);
	}
};

class IdentityFilter : public Filter {
public:
	IdentityFilter() {
		this->radius=0;
	}
	double f(double x) {
		return (x==0 ? 1 : 0);
	}
};
/*
class SobelHorizH : public Filter {
public:
	SobelHorizH() {
		this->radius=1;
	}
	double f(double x) {
		return (-x);
	}
};

class SobelHorizV : public Filter {
public:
	SobelHorizV() {
		this->radius=1;
	}
	double f(double x) {
		return (x==0 ? 2 : 1);
	}
};
*/


//convolve separable filter f of radius r with image I; basic algorithm from text

void ConvolveSep(R2Image *I, Filter *yfilter, Filter *xfilter) {
	int r_x = (int)xfilter->radius;
	int r_y = (int)yfilter->radius;
	int n_x = I->Width();
	int n_y = I->Height();
	R2Pixel *S = new R2Pixel[n_x];
	R2Image *I_out = new R2Image(n_x-2*r_x,n_y-2*r_y);
	//double d = I_out->Width();
	//I_out->Pixel(0x59,0x50).SetRed(3);
	//delete(I_out);
	for (int y=r_y; y<n_y-r_y; y++) {
		for (int x=0; x<n_x; x++) {
			S[x].Reset(0,0,0,0);
			//double d=filter->f(double(-r));
			//R2Pixel p =I->Pixel(x,y+r);
			for (int i=-r_y; i<=r_y; i++)  {
				//double d=filter->f(double(i));
				S[x] += yfilter->f(double(i))*I->Pixel(x,y-i);
				//				S[x].SetAlpha(filter->f(double(i))*I->Pixel(x,y-i).Alpha());
			}
		}
		for (int x=r_x; x<n_x-r_x; x++)
			for (int i=-r_x; i<=r_x; i++) {
				//printf("%x;%x\n",x,y);
				I_out->Pixel(x-r_x,y-r_y)+=S[x-i]*xfilter->f(double(i));
			}
	}
	for (int x=r_x; x<n_x-r_x; x++)
		for (int y=r_y; y<n_y-r_y; y++)
		{I->Pixel(x,y) = I_out->Pixel(x-r_x,y-r_y);}// printf("%x;%x\n",x-r,y-r);}
		//return I_out;
		delete(I_out);
}

void ConvolveSep(R2Image *I, Filter *filter) {
	ConvolveSep(I,filter,filter);
}

void ConvolveRowVec(R2Image *I, double *rowvec, int r) {
	int n_x = I->Width();
	int n_y = I->Height();
	R2Image *I_out = new R2Image(n_x-2*r,n_y);
	int c = r;
	double sum = 0;
	for (int i=0;i<2*r+1;i++) sum+=rowvec[i];
	if (sum==0) sum=1;
	for (int y=0; y<n_y; y++) 
		for (int x=r; x<n_x-r; x++) {
			for (int i=-r;i<=r;i++) 
				I_out->Pixel(x-r,y) += I->Pixel(x+i,y)*double(rowvec[c+i]);
			I_out->Pixel(x-r,y) /= sum;
		}



		for (int x=r; x<n_x-r; x++)
			for (int y=0; y<n_y; y++)
			{I->Pixel(x,y) = I_out->Pixel(x-r,y);}// printf("%x;%x\n",x-r,y-r);}
			//return I_out;
			delete(I_out);

}

void ConvolveColVec(R2Image *I, double *colvec, int r) {
	int n_x = I->Width();
	int n_y = I->Height();
	R2Image *I_out = new R2Image(n_x,n_y-2*r);
	int c = r;
	double sum = 0;
	for (int i=0;i<2*r+1;i++) sum+=colvec[i];
	if (sum==0) sum=1;
	for (int x=0; x<n_x; x++) 
		for (int y=r; y<n_y-r; y++) 
		{
			for (int j=-r;j<=r;j++) 
				I_out->Pixel(x,y-r) += I->Pixel(x,y+j)*double(colvec[c+j]);
			I_out->Pixel(x,y-r) /= sum;
		}

		for (int x=0; x<n_x; x++)
			for (int y=r; y<n_y-r; y++)
			{I->Pixel(x,y) = I_out->Pixel(x,y-r);}// printf("%x;%x\n",x-r,y-r);}
			//return I_out;
			delete(I_out);
}


//convolve separable filter f of radius r with image I; basic algorithm from text
/*void ConvolveSep(R2Image *I, Filter *filter) {
int r = filter->radius;
int n_x = I->Width();
int n_y = I->Height();
R2Pixel *S = new R2Pixel[n_x];
R2Image *I_out = new R2Image(n_x-2*r,n_y-2*r);
//double d = I_out->Width();
//I_out->Pixel(0x59,0x50).SetRed(3);
//delete(I_out);
for (int y=r; y<n_y-r; y++) {
for (int x=0; x<n_x; x++) {
S[x].Reset(0,0,0,0);
//double d=filter->f(double(-r));
//R2Pixel p =I->Pixel(x,y+r);
for (int i=-r; i<=r; i++)  {
//double d=filter->f(double(i));
S[x] += filter->f(double(i))*I->Pixel(x,y-i);
//				S[x].SetAlpha(filter->f(double(i))*I->Pixel(x,y-i).Alpha());
}
}
for (int x=r; x<n_x-r; x++)
for (int i=-r; i<=r; i++) {
//printf("%x;%x\n",x,y);
I_out->Pixel(x-r,y-r)+=S[x-i]*filter->f(double(i));
}
}
for (int x=r; x<n_x-r; x++)
for (int y=r; y<n_y-r; y++)
{I->Pixel(x,y) = I_out->Pixel(x-r,y-r);}// printf("%x;%x\n",x-r,y-r);}
//return I_out;
delete(I_out);
}*/


void R2Image::
Blur(double sigma)
{
	// Blur an image with a Gaussian filter with a given sigma.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Blur(%g) not implemented\n", sigma);
	GaussianFilter *g = new GaussianFilter(sigma);
	ConvolveSep(this,g);
}


void R2Image::
Sharpen()
{
	// Sharpen an image using a linear filter

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Sharpen() not implemented\n");
	double sigma=1;
	IdentityFilter d = IdentityFilter(); GaussianFilter g = GaussianFilter(sigma);
	double alpha=.1;
	SumFilter s = (d*(1+alpha)) + (g*(-alpha));
	//double k=(g*0).f(3);
	ConvolveSep(this,&s);
}

R2Image* ConvolveMat(R2Image *I,double *mat,int r) {
	R2Image *I_out = new R2Image(I->Width()-2*r,I->Height()-2*r);
	for (int x=r;x<I->Width()-r;x++) {
		for (int y=0;y<I->Height()-r;y++) {
			for (int i=-r;i<=r;i++) {
				for (int j=-r;j<=r;j++) {
					I_out->Pixel(x-r,y-r) += I->Pixel(x-i,y-j)*mat[(r+i)*(2*r+1)+(r+j)];
				}
			}
			I_out->Pixel(x-r,y-r).Clamp();
		}
	}
	return I_out;

	/*for (int x=r; x<I->Width()-r; x++)
	for (int y=r; y<I->Height()-r; y++)
	{I->Pixel(x,y) = I_out->Pixel(x-r,y-r);}// printf("%x;%x\n",x-r,y-r);}
	//return I_out;
	//delete(I_out);
	*/
}

void R2Image::
EdgeDetect(void)
{
	// Detect edges in an image.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "EdgeDetect() not implemented\n");
	//IdentityFilter d = IdentityFilter();
	//SobelHorizH shh = SobelHorizH();
	//SobelHorizV shv = SobelHorizV();
	//double shr[3] = {1,0,-1};
	//double shc[3] = {1,2,1};
	//ConvolveRowVec(this,shr,1);
	//ConvolveColVec(this,shc,1);
	double mat[] = {
		-1,-1,-1,
		-1,8,-1,
		-1,-1,-1
		};
	int r=1;
	//for (int i=-r;i<=r;i++) {
	//for (int j=-r;j<=r;j++) {
	//	printf("(%d,%d): %d\n",i,j,mat[(r-i)*(2*r+1)+(r-j)]);
	//}}
	*this=*(ConvolveMat(this,mat,r));
}



// Non-Linear filtering ////////////////////////////////////////////////


int compare (const void * a, const void * b)
{
	if (*(double*)a < *(double*)b) return -1;
	else
		if (*(double*)a > *(double*)b) return 1;
		else return 0;
}


void R2Image::
MedianFilter(double sigma)
{
	// Perform median filtering with a given width

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "MedianFilter(%g) not implemented\n", sigma);
	/*	double test[5] = {1,53.283,7.3,334,4};
	qsort(test,5,sizeof(double),compare);
	for (int i=0;i<5;i++) printf("%g\n",test[i]);
	*/
	/*	R2Image *I_out = new R2Image(this->width,this->height);

	double *nbhd_red;
	double *nbhd_green;
	double *nbhd_blue;
	int nbhd_size;

	r = int(sigma);

	for (int x=0; x<this->width; x++)
	for (int y=0; y<this->height; y++) {
	nbhd_size=0;
	nbhd_red = (double*)malloc(0);
	//nbhd_green = (double*)malloc(0);
	//nbhd_blue = (double*)malloc(0);
	for (int i=-min(x,sigma); i<=min(this->width-x-1,sigma); i++)
	for (int j=-min(y,sigma); j<=min(this->height-y-1,sigma); j++) {
	nbhd_size++;
	nbhd_red = (double*)realloc(nbhd_red,nbhd_size*sizeof(double));
	nbhd_red[nbhd_size-1] = this->Pixel(x+i,y+j).Red();
	//	printf("%d,%d; %d,%d\n",x,y,i,j);
	}
	//printf("%d,%d\n",x,y);
	//qsort(nbhd_red,nbhd_size,sizeof(double),compare);
	I_out->Pixel(x,y).SetRed(nbhd_red[int((nbhd_size+1)/2)-1]);
	free(nbhd_red);
	}
	*/
	R2Image *I_out = new R2Image(this->width,this->height);

	int r = int(sigma);

	int nbhd_size = (2*r+1)*(2*r+1);
	double *nbhd_red = new double[nbhd_size];
	double *nbhd_green  = new double[nbhd_size];
	double *nbhd_blue = new double[nbhd_size];


	for (int x=r; x<this->width-r; x++)
		for (int y=r; y<this->height-r; y++) {
			//nbhd_green = (double*)malloc(0);
			//nbhd_blue = (double*)malloc(0);
			for (int i=-r; i<=r; i++)
				for (int j=-r; j<=r; j++) {
					nbhd_red[(i+r)*(2*r+1)+j+r] = this->Pixel(x+i,y+j).Red();//printf("%g:",this->Pixel(x+i,y+j).Red());
					nbhd_green[(i+r)*(2*r+1)+j+r] = this->Pixel(x+i,y+j).Green();
					nbhd_blue[(i+r)*(2*r+1)+j+r] = this->Pixel(x+i,y+j).Blue();
				}
				//for (int d=0;d<nbhd_size;d++) printf("%g\n",nbhd_red[d]);printf("---\n");
				qsort(nbhd_red,nbhd_size,sizeof(double),compare);
				//for (int d=0;d<nbhd_size;d++) printf("%g\n",nbhd_red[d]); printf("---\n");
				qsort(nbhd_green,nbhd_size,sizeof(double),compare);
				qsort(nbhd_blue,nbhd_size,sizeof(double),compare);
				I_out->Pixel(x-r,y-r).SetRed(nbhd_red[int((nbhd_size)/2)-1]);
				I_out->Pixel(x-r,y-r).SetGreen(nbhd_green[int((nbhd_size)/2)-1]);
				I_out->Pixel(x-r,y-r).SetBlue(nbhd_blue[int((nbhd_size)/2)-1]);
				//free(nbhd_red);
		}

		for (int x=r; x<this->width-r; x++)
			for (int y=r; y<this->height-r; y++)
				this->Pixel(x,y) = I_out->Pixel(x-r,y-r);// printf("%x;%x\n",x-r,y-r);}
		//return I_out;
		delete(I_out);
}

void R2Image::
BilateralFilter(double rangesigma, double domainsigma)
{
	// Perform bilateral filtering with a given range and domain widths.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "BilateralFilter(%g, %g) not implemented\n", rangesigma, domainsigma);
}


// Resampling operations  ////////////////////////////////////////////////

class BoxFilter : public Filter {
public:
	BoxFilter(double radius) {
		this->radius = radius;
	}
	double f(double x) {
		return (abs(x)>radius ? 0 : (double(1)/(2*radius)));
	}
	double f(double x,double y) {
		return this->f(x) * this->f(y);
	}
};

class TentFilter : public Filter {
public:
	TentFilter(int radius) {
		this->radius = radius;
	}
	double f(double x) {
		return (abs(x)>radius ? 0 : double(radius-abs(x))/radius);
	}
	double f(double x,double y) {
		return this->f(x) * this->f(y);
	}

};


R2Pixel reconstruct_x (R2Image *I, Filter *filter, double u, int y) {
	R2Pixel p = R2Pixel();
	for (int i=(int)max(0,ceil(u-filter->radius)); i<=(int)min(I->Width()-1,floor(u+filter->radius)); i++) {
		p += I->Pixel(i,y)*filter->f(u-double(i));
		//printf("%d;%d\n",i,y);
	}
	return p;
}

R2Pixel reconstruct_y (R2Image *I, Filter *filter, int x, double v) {
	R2Pixel p = R2Pixel();
	for (int j=(int)max(0,ceil(v-filter->radius)); j<=(int)min(I->Height()-1,floor(v+filter->radius)); j++) {
		p += I->Pixel(x,j)*filter->f(v-double(j));
		//printf("%x;%x\n",x,j);
	}
	return p;
}

R2Pixel reconstruct (R2Image *I, Filter *filter, double u, double v) {
	R2Pixel p = R2Pixel(0,0,0,1);

	for (int i=(int)max(0,ceil(u-filter->radius)); i<=(int)min(I->Width()-1,floor(u+filter->radius)); i++) {
	for (int j=(int)max(0,ceil(v-filter->radius)); j<=(int)min(I->Height()-1,floor(v+filter->radius)); j++) {
		p += I->Pixel(i,j)*filter->f(u-double(i),v-double(j));
		//printf("%x;%x\n",x,j);
	}
	}
	return p;	
}

R2Pixel reconstruct (R2Image *I, Filter *filter, R2Point p0) {
	return reconstruct(I,filter,p0.X(),p0.Y());
}


void R2Image::
Scale(double sx, double sy, int sampling_method)
{
	// Scale an image in x by sx, and y by sy.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Scale(%g, %g, %d) not implemented\n", sx, sy, sampling_method);

	Filter *filter;

	if (sampling_method==0) filter = new BoxFilter(.5);
	else if (sampling_method==1)  filter = new TentFilter(1);
	else filter = new GaussianFilter(.5*sqrt(sx*sy));


	int nx_out = (int)floor(sx*this->width); 
	int ny_out = (int)floor(sy*this->height);
	double u,v;

	R2Image *I_out = new R2Image(nx_out,ny_out);//ff,117; width==0x55=85
	for (int j=0;j<ny_out;j++) {
		for (int i=0;i<nx_out;i++) {
			u = double(i/sx);
			v = double(j/sy);
			I_out->Pixel(i,j) = reconstruct(this,filter,u,v);
			//printf("%x,%x\n",i,j);
		}
	}

	*this = *I_out;
	//below is an implementation scaling first the x direction, then the y direction

/*
	R2Image *I_out_x = new R2Image(nx_out,this->height);
	for (int j=0;j<this->height;j++) {
		for (int i=0;i<nx_out;i++) {
			u = double(i/sx);
			I_out_x->Pixel(i,j) = reconstruct_x(this,filter,u,j);
		}
	}

	R2Image *I_out_y = new R2Image(nx_out,ny_out);

	for (int i=0;i<nx_out;i++) {
		for (int j=0;j<ny_out;j++) {
			v = double(j/sy);
			I_out_y->Pixel(i,j) = reconstruct_y(I_out_x,filter,i,v);
		}
	}
	delete(I_out_x);
	*this = *I_out_y;
*/
}

R2Point *rotate(R2Point *p, double angle) {
	return new R2Point(cos(angle)*p->X() - sin(angle)*p->Y(),
						sin(angle)*p->X() + cos(angle)*p->Y());
}

#define sqr(a) ((a)*(a))
double distance(R2Point *a,R2Point *b) {
	return sqrt(sqr(a->X()-b->X())+sqr(a->Y()-b->Y()));
}

void R2Image::
Rotate(double angle, int sampling_method)
{
	// Rotate an image by the given angle.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Rotate(%g, %d) not implemented\n", angle, sampling_method);


	Filter *filter;
	if (sampling_method==0) filter = new BoxFilter(.5);
	else if (sampling_method==1)  filter = new TentFilter(1);
	else filter = new GaussianFilter(.5);

	double pi = 3.1416;
	angle = (angle/abs(angle))*fmod(abs(angle),2.0*pi);

	R2Point *topleft = rotate(new R2Point(0,this->height-1),angle);
	R2Point *topright = rotate(new R2Point(this->width-1,this->height-1),angle);
	R2Point *bottomleft = rotate(new R2Point(0,0),angle);
	R2Point *bottomright = rotate(new R2Point(this->width-1,0),angle);
	
#define min4(a,b,c,d) min(a,min(b,min(c,d)))
#define max4(a,b,c,d) max(a,max(b,max(c,d)))
	double leftedge = min4(topleft->X(),topright->X(),bottomleft->X(),bottomright->X());
	double rightedge = max4(topleft->X(),topright->X(),bottomleft->X(),bottomright->X());
	double topedge = max4(topleft->Y(),topright->Y(),bottomleft->Y(),bottomright->Y());
	double bottomedge = min4(topleft->Y(),topright->Y(),bottomleft->Y(),bottomright->Y());

	int ny = (int)ceil(topedge-bottomedge+1);
	int nx = (int)ceil(rightedge-leftedge+1);
	R2Image *I_out = new R2Image(nx,ny);

	R2Point *p;
	for (int x=0; x<nx; x++) {
		for (int y=0; y<ny; y++) {
			p = rotate(new R2Point(x+leftedge,y+bottomedge),-angle);
			if ((p->X() < this->width) && (p->Y() < this->height))
				I_out->Pixel(x,y) = reconstruct(this,filter,p->X(),p->Y());
		}
	}

	*this = *I_out;
}


void R2Image::
MotionBlur(R2Vector& movement, int sampling_method)
{
	// Perform Motion Blur as if the camera translated along the given vector

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//only looks good for motion parallel to axes, per the assignment
	//fprintf(stderr, "MotionBlur(%g %g) not implemented\n", movement.X(), movement.Y(), sampling_method);
	int rx=int(movement.X());
	int ry=int(movement.Y());
	if (rx==0) rx=1; if (ry==0) ry=1;
	double *mbx = new double[2*rx+1];
	double *mby = new double[2*ry+1];
	for (int i=0;i<rx;i++) mbx[i]=0;
	for (int i=0;i<ry;i++) mby[i]=0;
	for (int i=rx;i<2*rx+1;i++) mbx[i]=1;
	for (int i=ry;i<2*ry+1;i++) mby[i]=1;
	ConvolveRowVec(this,mbx,rx);
}



void R2Image::
Fun(int sampling_method)
{
	// Warp an image using a creative filter of your choice.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Fun(%d) not implemented\n", sampling_method);
	
	Filter *filter;
	if (sampling_method==0) filter = new BoxFilter(.5);
	else if (sampling_method==1)  filter = new TentFilter(1);
	else filter = new GaussianFilter(1);

	double angle = 2;

	R2Point center = R2Point((this->width-1)/2,(this->height-1)/2);
	R2Image *I_out = new R2Image(this->width,this->height);
//	R2Point *p0;
	double u,v;

	for (int x=0; x<this->width; x++)
		for (int y=0; y<this->height; y++) {
			//p0 = rotate(new R2Point(x,y),angle/distance(new R2Point(x,y),&center));
			//I_out->Pixel(x,y) = reconstruct(this,filter,*p0);
			u = abs(x-center.X())*angle;
			v = abs(y-center.Y())*angle;
			I_out->Pixel(x,y) = reconstruct(this,filter,u,v);
		}

		*this = *I_out;
}


// Dither operations ////////////////////////////////////////////////

void R2Image::
Quantize (int nbits)
{
	// Quantizes an image with "nbits" bits per channel.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Quantize(%d) not implemented\n", nbits);
	double nbuckets = pow(2.0,nbits);
	double delta = double(nbuckets)/10000;
	for (int x=0; x<this->width; x++)
		for (int y=0; y<this->height; y++) {
			this->Pixel(x,y).SetRed(floor(this->Pixel(x,y).Red()*(double(nbuckets)-delta))/(double(nbuckets)-delta));
			this->Pixel(x,y).SetGreen(floor(this->Pixel(x,y).Green()*(double(nbuckets)-delta))/(double(nbuckets)-delta));
			this->Pixel(x,y).SetBlue(floor(this->Pixel(x,y).Blue()*(double(nbuckets)-delta))/(double(nbuckets)-delta));
			this->Pixel(x,y) *= double(nbuckets)/(nbuckets-1); //scale to fill out 0..1
		}

}



void R2Image::
RandomDither(int nbits)
{
	// Converts and image to nbits per channel using random dither.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "RandomDither(%d) not implemented\n", nbits);
	//first run through image adding noise proportional to bucket size
		double nbuckets = pow(2.0,nbits);

		for (int x=0; x<width; x++)
		for (int y=0; y<height; y++)
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS-1; k++){
				double sign = (rand() < RAND_MAX/2 ? 1 : 1);
				Pixel(x,y).Components()[k]+=sign*double(rand())/RAND_MAX/nbuckets;
					Pixel(x,y).Clamp();
			}

			//then quantize
			this->Quantize(nbits);
}



void R2Image::
OrderedDither(int nbits)
{
	// Converts an image to nbits per channel using ordered dither, 
	// with a 4x4 Bayer's pattern matrix.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "OrderedDither(%d) not implemented\n", nbits);
	int Bayer[2][2] = {1,3,
						4,2};
	int i,j;
	//algorithm, notation from slides
	for (int x=0; x<this->width; x++)
		for (int y=0; y<this->height; y++) {
				i = x % 2;
				j = y % 2;
				this->Pixel(x,y) *= pow(2.0,nbits)-1;
				for (int k=0; k<R2_IMAGE_NUM_CHANNELS; k++) {
					double c = this->Pixel(x,y).Components()[k];
					if (c-floor(c) < double(Bayer[i][j])/5.0)
						Pixel(x,y).Components()[k] = floor(c)/(pow(2.0,nbits)-1);
					else
						Pixel(x,y).Components()[k] = ceil(c)/(pow(2.0,nbits)-1);
				}
				
		}


}



void R2Image::
FloydSteinbergDither(int nbits)
{
	// Converts an image to nbits per channel using Floyd-Steinberg dither.
	// with error diffusion.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "FloydSteinbergDither(%d) not implemented\n", nbits);
	double nbuckets = pow(2.0,nbits);
	double oldval, error;
	double delta = 0;//double(nbuckets)/10000;
	for (int y=this->height-1; y>0; y--)
		for (int x=1; x<this->width-1; x++) {
			for (int k=0; k<R2_IMAGE_NUM_CHANNELS; k++) {
				oldval = this->Pixel(x,y).Components()[k];
				double c;
				this->Pixel(x,y).Components()[k] = (c=floor(oldval*(double(nbuckets)-delta))/(double(nbuckets)-delta));
				this->Pixel(x,y).Components()[k] *= double(nbuckets)/(nbuckets-1);

				error = oldval - Pixel(x,y).Components()[k];
				this->Pixel(x+1,y).Components()[k] += (7.0/16.0)*error;
				this->Pixel(x+1,y-1).Components()[k] += (1.0/16.0)*error;
				this->Pixel(x,y-1).Components()[k] += (5.0/16.0)*error;
				this->Pixel(x-1,y-1).Components()[k] += (3.0/16.0)*error;
			}
		}
}



// Miscellaneous operations ////////////////////////////////////////////////

void R2Image::
Crop(int x, int y, int w, int h)
{
	// Extracts a sub image from the image, 
	// at position (x, y), width w, and height h.

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Crop(%d %d  %d %d) not implemented\n", x, y, w, h);

	assert(x>=0 && y>=0); 
	assert(x+w<=this->width); 
	assert(y+h<=this->height);

	R2Image *I_out = new R2Image(w,h);
	for (int i=0; i<w; i++)
		for (int j=0; j<h; j++)
			I_out->Pixel(i,j) = this->Pixel(x+i,y+j);

	*this = *I_out;
}

void R2Image::
Composite(const R2Image& top, int operation)
{
	// Composite passed image on top of this one using operation (e.g., OVER)

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "Composite not implemented\n");
}

void R2Image::
ExtractChannel(int channel)
{
	// Extracts a channel of an image (e.g., R2_IMAGE_RED_CHANNEL).  
	// Leaves the specified channel intact, 
	// and sets all the other ones to zero.

	// MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
	//fprintf(stderr, "ExtractChannel(%d) not implemented\n", channel);
	for (int k=0; k<R2_IMAGE_NUM_CHANNELS; k++) 
		if (k!=channel)
			for (int i=0; i<width; i++)
				for (int j=0; j<height; j++)
					Pixel(i,j).Components()[k] = 0;
}


void R2Image::
CopyChannel(const R2Image& from_image, int from_channel, int to_channel)
{
	// Copies one channel of an image (e.g., R2_IMAGE_RED_CHANNEL).  
	// to another channel

	// MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
	fprintf(stderr, "CopyChannel not implemented\n");
}

void R2Image::
Add(const R2Image& image)
{
	// Add passed image pixel-by-pixel.

	// MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
	fprintf(stderr, "Add not implemented\n");
}

void R2Image::
Subtract(const R2Image& image)
{
	// Subtract passed image from this image.

	// MAY FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE) (NO CREDIT FOR ASSIGNMENT)
	fprintf(stderr, "Subtract not implemented\n");
}

//a few helper functions for morphing...

R2Image *dissolve(R2Image *src, R2Image *dest, double t) {
	int nx = min(src->Width(),dest->Width());
	int ny = min(src->Height(),dest->Height());
	R2Image *I_out = new R2Image(nx,ny);
	
	for (int x=0; x<nx; x++)
	for (int y=0; y<ny; y++) {
			I_out->Pixel(x,y) = (1-t)*src->Pixel(x,y) + t*dest->Pixel(x,y);
	}

	return I_out;
}

R2Point morphPixel(R2Image *src, R2Segment source_segment, R2Segment target_segment, R2Point &X) {
	R2Vector PX = X.Vector()-R2Vector(target_segment.Start().X(),target_segment.Start().Y());
	//u==(X-P).(Q-P)/length(Q-P)
	double u = PX.Dot(target_segment.Vector());
	u /= target_segment.Length();
	//v==(X-P).Perpendicular(Q-P)/length(Q-P)
	double v = PX.Dot(target_segment.Normal());

	//X' == P' + u*(Q'-P') + v*Perp(Q'-P')/dist(Q'-P')
	R2Point X0 = source_segment.Start() + u*(source_segment.Vector()*source_segment.Length());
	X0 += v*source_segment.Normal();

	return X0;
}

double distPointSeg(R2Point &X,R2Segment &s) {
	R2Vector PX = X.Vector()-s.Start().Vector();
	R2Vector PQ = s.Vector(); 
	double u = PX.Dot(PQ);
	if (u<0)
		return PX.Length();
	else
		if (u>PQ.Length())
			return (X.Vector()-s.End().Vector()).Length();
		else {
			R2Point Xproj = X; Xproj.Project(s.Line());
			return (Xproj.Vector()-X.Vector()).Length();
		}
}

double morphWeight(double length,double dist,double a, double b, double p) {
	return pow((pow(length,p)/(a+dist)),b);
}

double morphWeight(double length,double dist) {
	return morphWeight(length,dist,1,1,1);
}

R2Pixel morphPixelMulti(R2Image *src, R2Segment *source_segment, R2Segment *target_segment, int nsegments,
				   Filter *filter, R2Point &X) {
		double dist,weight,weightsum=0;
		R2Vector D_i; //displacements X -> X_i'
		R2Point X0_i;
		R2Vector dsum = R2Vector(0,0);

		for (int i=0;i<nsegments;i++) {
			X0_i = morphPixel(src,source_segment[i],target_segment[i],X);
			D_i = X0_i.Vector() - X.Vector();
			dist = distPointSeg(X,target_segment[i]);
			weight = morphWeight(target_segment[i].Length(),dist);
			weightsum += weight;
			dsum += weight * D_i;
		}

		R2Point X0 = X+(dsum/weightsum);
		return reconstruct(src,filter,X0.X(),X0.Y());
}

R2Segment *interpolateSeg(R2Segment *seg1, R2Segment *seg2, double t) {
	return new R2Segment((1-t)*seg1->Start()+t*seg2->Start(),(1-t)*seg1->End()+t*seg2->End());
}
//mark segment in an image, for debugging purposes
void markSeg(R2Image *I,R2Segment seg,R2Pixel colorPixel) {
	int r = 3;
	for (int i=-r;i<=r;i++)
		for (int j=-r;j<=r;j++) {
			I->Pixel((int)seg.Start().X()+i,(int)seg.Start().Y()+j) = colorPixel;
			I->Pixel((int)seg.End().X()+i,(int)seg.End().Y()+j) = colorPixel;
		}
}

void R2Image::
Morph(const R2Image& target, 
	  R2Segment *source_segments, R2Segment *target_segments, int nsegments, 
	  double t, int sampling_method)
{
	// Morph this source image towards a passed target image by t using pairwise line segment correspondences

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	//fprintf(stderr, "Morph not implemented\n");


	Filter *filter;
	if (sampling_method==0) filter = new BoxFilter(.5);
	else if (sampling_method==1)  filter = new TentFilter(1);
	else filter = new GaussianFilter(1);

	R2Image *dest = new R2Image(target);
	//*this = *dissolve(this,dest,.5);

	//debugging lines
	/*for (int i=0;i<nsegments;i++) {
	markSeg(this,source_segments[i],R2green_pixel);
	markSeg(dest,target_segments[i],R2cyan_pixel);
	};//*this=*dest;
	//return;
	//double a = distPointSeg(R2Point(0,10),R2Segment(20,0,20,100));
*/

	R2Image *I0_out = new R2Image(this->Width(),this->Height());
	R2Image *I1_out = new R2Image(dest->Width(),dest->Height());


	//get interpolated line segments
	R2Segment *mid_segments = new R2Segment[nsegments];
	for (int i=0;i<nsegments;i++) mid_segments[i] = *interpolateSeg(&source_segments[i],&target_segments[i],t);

	//int i=0;

	//distort source image toward intermediate segments
	for (int x=0; x<I0_out->Width(); x++)
		for (int y=0; y<I0_out->Height(); y++) {
//			I0_out->Pixel(x,y) = morphPixel(this,source_segments,mid_segments,filter,R2Point(x,y));//(dest,filter,X0.X(),X0.Y());
			I0_out->Pixel(x,y) = morphPixelMulti(this,source_segments,mid_segments,nsegments,filter,R2Point(x,y));//(dest,filter,X0.X(),X0.Y());
		}

	//distort destination image toward intermediate segments
	for (int x=0; x<I1_out->Width(); x++)
		for (int y=0; y<I1_out->Height(); y++) {
			//I1_out->Pixel(x,y) = morphPixel(dest,target_segments,mid_segments,filter,R2Point(x,y));//(dest,filter,X0.X(),X0.Y());
			I1_out->Pixel(x,y) = morphPixelMulti(dest,target_segments,mid_segments,nsegments,filter,R2Point(x,y));//(dest,filter,X0.X(),X0.Y());

		}


	//cross dissolve two intermediate segments
		*this = *dissolve(I0_out,I1_out,t);

}


////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
	// Initialize everything
	if (pixels) { delete [] pixels; pixels = NULL; }
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
	unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
	unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
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
	assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

	// Assign width, height, and number of pixels
	width = bmih.biWidth;
	height = bmih.biHeight;
	npixels = width * height;

	// Allocate unsigned char buffer for reading pixels
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = bmih.biSizeImage;
	unsigned char *buffer = new unsigned char [nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Read buffer 
	fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
	if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
		fprintf(stderr, "Error while reading BMP file %s", filename);
		return 0;
	}

	// Close file
	fclose(fp);

	// Allocate pixels for image
	pixels = new R2Pixel [ width * height ];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double b = (double) *(p++) / 255;
			double g = (double) *(p++) / 255;
			double r = (double) *(p++) / 255;
			R2Pixel pixel(r, g, b, 1);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete [] buffer;

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
	bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
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
			fputc((unsigned char) b, fp);
			fputc((unsigned char) g, fp);
			fputc((unsigned char) r, fp);
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
	pixels = new R2Pixel [ width * height ];
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
		for (int j = height-1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				double r = (double) getc(fp) / max_value;
				double g = (double) getc(fp) / max_value;
				double b = (double) getc(fp) / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}
	else {
		// Read asci image data 
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height-1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				// Read pixel values
				int red, green, blue;
				if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
					fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
					fclose(fp);
					return 0;
				}

				// Assign pixel values
				double r = (double) red / max_value;
				double g = (double) green / max_value;
				double b = (double) blue / max_value;
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
		for (int j = height-1; j >= 0 ; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int) (255 * p.Red());
				int g = (int) (255 * p.Green());
				int b = (int) (255 * p.Blue());
				fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
				if (((i+1) % 4) == 0) fprintf(fp, "\n");
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
		for (int j = height-1; j >= 0 ; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int) (255 * p.Red());
				int g = (int) (255 * p.Green());
				int b = (int) (255 * p.Blue());
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
	pixels = new R2Pixel [ npixels ];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Allocate unsigned char buffer for reading image
	int rowsize = ncomponents * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char [nbytes];
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
				r = g = b = (double) *(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 1) {
				r = g = b = (double) *(p++) / 255;
				a = 1;
				p++;
			}
			else if (ncomponents == 3) {
				r = (double) *(p++) / 255;
				g = (double) *(p++) / 255;
				b = (double) *(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 4) {
				r = (double) *(p++) / 255;
				g = (double) *(p++) / 255;
				b = (double) *(p++) / 255;
				a = (double) *(p++) / 255;
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
	delete [] buffer;

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
	jpeg_set_quality(&cinfo, 75, TRUE);
	jpeg_start_compress(&cinfo, TRUE);

	// Allocate unsigned char buffer for reading image
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char [nbytes];
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
			int r = (int) (255 * pixel.Red());
			int g = (int) (255 * pixel.Green());
			int b = (int) (255 * pixel.Blue());
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
	delete [] buffer;

	// Return number of bytes written
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}






