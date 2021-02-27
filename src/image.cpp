//CSCI 5607 HW 2 - Image Conversion Instructor: S. J. Guy <sjguy@umn.edu>
//In this assignment you will load and convert between various image formats.
//Additionally, you will manipulate the stored image data by quantizing, cropping, and supressing channels

#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <fstream>
using namespace std;


float clamp(float n, float lower, float upper) {
  return std::max(lower, std::min(n, upper));
}

/**
 * Image
 **/
Image::Image (int width_, int height_){
    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
		int b = 0; //which byte to write to
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				data.raw[b++] = 0;
				data.raw[b++] = 0;
				data.raw[b++] = 0;
				data.raw[b++] = 255;
			}
		}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){
	width           = src.width;
	height          = src.height;
	num_pixels      = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;
	
	data.raw = new uint8_t[num_pixels*4];
	
	memcpy(data.raw, src.data.raw, sizeof(Pixel) * num_pixels);
	//*data.raw = *src.data.raw;
}

Image::Image (char* fname){

	int lastc = strlen(fname);
   int numComponents; //(e.g., Y, YA, RGB, or RGBA)
   data.raw = stbi_load(fname, &width, &height, &numComponents, 4);
	
	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}
	
	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;
	
}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){
	int lastc = strlen(fname);

	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}


void Image::Brighten (double factor){
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			Pixel scaled_p = p*factor;
			GetPixel(x,y) = scaled_p;
		}
	}
}

// Keep only non-zero red, green, or blue components for the channel value 0, 1, and 2 respectively
void Image::ExtractChannel(int channel) {
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			Pixel new_p;
			switch(channel) {
				case 1: {
					new_p = Pixel(0, p.g);
					break;
				}
				case 2: {
					new_p = Pixel(0, 0, p.b);
					break;
				}
				case 3: {
					new_p = Pixel(0, 0, 0, p.a);
					break;
				}
				case 0:
				default: {
					new_p = Pixel(p.r);
					break;
				}
			}
			GetPixel(x,y) = new_p;
		}
	}
}


// Quantize the intensities stored for each pixel's values into 2^nbits possible equally-spaced values
void Image::Quantize (int nbits) {
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			GetPixel(x,y) = PixelQuant(p, nbits);
		}
	}
}

//Crop an image to a rectangle starting at (x,y) with a width w and a height h
Image* Image::Crop(int x0, int y0, int w, int h) {
	Image* new_image = new Image(w, h);
	for (int x = 0; x < min(w, x0 + Width()); x++){
		for (int y = 0; y < min(h, y0 + Height()); y++){
			Pixel p = GetPixel(x0 + x, y0 + y);
			new_image->GetPixel(x,y) = p;
		}
	}
	return new_image;
}


void Image::AddNoise (double factor){
	// For each pixel, create a random pixel, then lerp to it by the factor.
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			GetPixel(x,y) = PixelLerp(p, PixelRandom(), factor);
		}
	}
}

void Image::ChangeContrast (double factor){
	double luminance = 0;
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			luminance += p.Luminance();
		}
	}

	luminance /= Width() * Height();

	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			Pixel average_l = PixelLerp(Pixel(0,0,0), Pixel(255,255,255), luminance);
			GetPixel(x, y) = PixelLerp(p, average_l, 1-factor);
		}
	}
}


void Image::ChangeSaturation(double factor){
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			GetPixel(x, y) = PixelLerp(Grayscale(255.0 * p.Luminance()), p, factor);
		}
	}
}

//For full credit, check that your dithers aren't making the pictures systematically brighter or darker
// Adding linear noise with a max of the quantize delta would achieve the same, I think.
void Image::RandomDither (int nbits){
	int x,y;
	for (x = 0 ; x < Width() ; x++){
		for (y = 0 ; y < Height() ; y++){
			Pixel p = GetPixel(x, y);
			int shift = 8-nbits;
			float mult = 255.0/(255 >> shift);

			int new_r = (p.r >> shift);
			int new_g = (p.g >> shift);
			int new_b = (p.b >> shift);

			float ran = (((float) rand()) / INT_MAX);

			float rd = ((float) (new_r * mult) - p.r) / ((0b1 << shift) - 1);
			float gd = ((float) (new_g * mult) - p.g) / ((0b1 << shift) - 1);
			float bd = ((float) (new_b * mult) - p.b) / ((0b1 << shift) - 1);

			if ((rd < 0) && (-2*ran > rd)) new_r++;
			if ((gd < 0) && (-2*ran > gd)) new_g++;
			if ((bd < 0) && (-2*ran > bd)) new_b++;

			if ((rd > 0) && (2*ran < rd)) new_r--;
			if ((gd > 0) && (2*ran < gd)) new_g--;
			if ((bd > 0) && (2*ran < bd)) new_b--;

			GetPixel(x, y).SetClamp(new_r*mult , new_g*mult, new_b*mult);
		}
	}
}

//This bayer method gives the quantization thresholds for an ordered dither.
//This is a 4x4 dither pattern, assumes the values are quantized to 16 levels.
//You can either expand this to a larger bayer pattern. Or (more likely), scale
//the threshold based on the target quantization levels.
static int Bayer4[4][4] ={
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits){
	/* WORK HERE */
}

/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits){
	for (int x = 0 ; x < Width(); x++) {
		for (int y = 0 ; y < Height(); y++) {
			Pixel p = GetPixel(x, y);
			int shift = 8-nbits;
			float mult = 255.0/(255 >> shift);

			int qr = mult*(p.r >> shift);
			int qg = mult*(p.g >> shift);
			int qb = mult*(p.b >> shift);

			int rd = p.r - qr;
			int gd = p.g - qg;
			int bd = p.b - qb;

			if (ValidCoord(x+1, y)) {
				Pixel np = GetPixel(x+1, y);
				np.SetClamp(np.r + (rd * ALPHA), np.g + (gd * ALPHA), np.b + (bd * ALPHA));
				GetPixel(x+1, y) = np;
			}
			if (ValidCoord(x-1, y+1)) {
				Pixel np = GetPixel(x-1, y+1);
				np.SetClamp(np.r + (rd * ALPHA), np.g + (gd * ALPHA), np.b + (bd * ALPHA));
				GetPixel(x-1, y+1) = np;
			}
			if (ValidCoord(x, y+1)) {
				Pixel np = GetPixel(x, y+1);
				np.SetClamp(np.r + (rd * ALPHA), np.g + (gd * ALPHA), np.b + (bd * ALPHA));
				GetPixel(x, y+1) = np;
			}
			if (ValidCoord(x+1, y+1)) {
				Pixel np = GetPixel(x+1, y+1);
				np.SetClamp(np.r + (rd * ALPHA), np.g + (gd * ALPHA), np.b + (bd * ALPHA));
				GetPixel(x+1, y+1) = np;
			}
			GetPixel(x, y).SetClamp(qr, qg, qb);
		}
	}
}

void Image::Blur(int n){
	if (n <= 0) { return; }

	vector<float> filter(n+1);

	float coef = -1/n;
	for (int i = 0; i < filter.size(); i++) {
		filter[i] = exp(coef*pow(i, 2));
	}
	float sum = filter[0];
	for (int i = 1; i < filter.size(); i++) { sum += 2*filter[i]; }
	for (int i = 0; i < filter.size(); i++) { filter[i] /= sum; }

	for (int dim = 0; dim <= 1; dim++) {
		Image* img_copy = new Image(*this);
		for (int x = 0 ; x < Width() ; x++){
			for (int y = 0 ; y < Height() ; y++){
				float rt=0, gt=0, bt=0;
				for (int i = -n; i <= n; i++) {
					Pixel p = img_copy->GetPixelClamp(x+i*(1-dim), y+i*dim);
					rt += filter[abs(i)] * pow(p.r,2);
					gt += filter[abs(i)] * pow(p.g,2);
					bt += filter[abs(i)] * pow(p.b,2);
				}
				GetPixel(x, y).SetClamp(sqrt(rt), sqrt(gt), sqrt(bt));
			}
		}
		delete img_copy;
	}
}

void Image::Sharpen(int n){
	Image* blurred = new Image(*this);
	blurred->Blur(n);
	for (int x = 0; x < Width(); x++) {
		for (int y = 0; y < Height(); y++) {
			Pixel bp = blurred->GetPixel(x, y);
			Pixel p = GetPixel(x, y);
			GetPixel(x,y) = PixelLerp(bp, p, 2.0);
		}
	}
}

static int gm[3][3] ={
    {1, 0, -1},
    {2, 0, -2},
    {1, 0, -1},
};

void Image::EdgeDetect(){
	Image gx(Width(), Height());
	Image gy(Width(), Height());

	// detect x
	for (int x = 0; x < Width(); x++) {
		for (int y = 0; y < Height(); y++) {
			int rr=0, gg=0, bb=0;
			for (int xi = -1; xi <= 1; xi++) {
				for (int yi = -1; yi <= 1; yi++) {
					rr += GetPixelClamp(x+xi, y+yi).r * gm[yi+1][xi+1];
					gg += GetPixelClamp(x+xi, y+yi).g * gm[yi+1][xi+1];
					bb += GetPixelClamp(x+xi, y+yi).b * gm[yi+1][xi+1];
				}
			}
			gx.GetPixel(x, y).SetClamp(abs(rr), abs(gg), abs(bb));
		}
	}
	// detect y
	for (int x = 0; x < Width(); x++) {
		for (int y = 0; y < Height(); y++) {
			int rr=0, gg=0, bb=0;
			for (int xi = -1; xi <= 1; xi++) {
				for (int yi = -1; yi <= 1; yi++) {
					rr += GetPixelClamp(x+xi, y+yi).r * gm[xi+1][yi+1];
					gg += GetPixelClamp(x+xi, y+yi).g * gm[xi+1][yi+1];
					bb += GetPixelClamp(x+xi, y+yi).b * gm[xi+1][yi+1];
				}
			}
			gy.GetPixel(x, y).SetClamp(abs(rr), abs(gg), abs(bb));
		}
	}
	// combine
	for (int x = 0; x < Width(); x++) {
		for (int y = 0; y < Height(); y++) {
			Pixel px = gx.GetPixel(x, y);
			Pixel py = gy.GetPixel(x, y);
			int red = sqrt(pow(px.r, 2) + pow(py.r, 2));
			int green = sqrt(pow(px.g, 2) + pow(py.g, 2));
			int blue = sqrt(pow(px.b, 2) + pow(py.b, 2));
			GetPixel(x,y).SetClamp(red, green, blue);
		}
	}
}

Image* Image::Scale(double sx, double sy){
	// sampling method
	Image* out_image = new Image(Width() * sx, Height() * sy);
	
	for (int out_x = 0; out_x < out_image->width; out_x++) {
		for (int out_y = 0; out_y < out_image->height; out_y++) {
			if (sampling_method == IMAGE_SAMPLING_POINT) {
				out_image->GetPixel(out_x, out_y) = GetPixel((int) out_x / sx, (int) out_y / sy);
			}
			if (sampling_method == IMAGE_SAMPLING_BILINEAR) {
				if (sx < 1.0 && sy < 1.0) {
					float sc_x = out_x / sx;
					float sc_y = out_y / sy;

					float sum = 0;
					int r=0, g=0, b=0;
					for (int samx = ceil(0.001+ sc_x - 1/sx); samx < (int) sc_x + 1/sx; samx++) {
						for (int samy = ceil(0.001+ sc_y - 1/sy); samy < (int) sc_y + 1/sy; samy++) {
							float factor = 1/sx + 1/sy - abs(sc_x - samx) - abs(sc_y - samy);
							sum += factor;

							Pixel p = GetPixelClamp(samx, samy);
							r += p.r * factor;
							g += p.g * factor;
							b += p.b * factor;
						}
					}
					out_image->GetPixel(out_x, out_y) = Pixel(r/sum, g/sum, b/sum);
				} else if (sx > 1.0 && sy > 1.0) {
					float sc_x = out_x / sx;
					float sc_y = out_y / sy;
					Pixel p = GetPixelClamp((int) sc_x, (int) sc_y);
					Pixel px = GetPixelClamp((int) ceil(sc_x), (int) sc_y);
					Pixel py = GetPixelClamp((int) sc_x, (int) ceil(sc_y));
					Pixel pxy = GetPixelClamp((int) ceil(sc_x), (int) ceil(sc_y));

					Pixel pt = PixelLerp(p, px, 1-(ceil(sc_x) - sc_x));
					Pixel pb = PixelLerp(py, pxy, 1-(ceil(sc_x) - sc_x));
					out_image->GetPixel(out_x, out_y) = PixelLerp(pt, pb, 1-(ceil(sc_y) - sc_y));
				} else {
					printf("uneven isn't supported lol\n");
				}
			}
			if (sampling_method == IMAGE_SAMPLING_GAUSSIAN) {
				if (sx < 1.0 && sy < 1.0) {
					float sc_x = out_x / sx;
					float sc_y = out_y / sy;

					float sum = 0;
					int r=0, g=0, b=0;

					for (int samx = ceil(sc_x - 1/sx); samx <= (int) sc_x + 1/sx; samx++) {
						for (int samy = ceil(sc_y - 1/sy); samy <= (int) sc_y + 1/sy; samy++) {
							float x_dist = abs(samx - sc_x);
							float y_dist = abs(samy - sc_y);

							float factor = pow(x_dist, 2) + pow(y_dist, 2);
							factor = exp(-0.5 * sx * factor);
							sum += factor;

							Pixel p = GetPixelClamp(samx, samy);
							r += p.r * factor;
							g += p.g * factor;
							b += p.b * factor;
						}
					}
					out_image->GetPixel(out_x, out_y) = Pixel(r/sum, g/sum, b/sum);
				} else if (sx > 1.0 && sy > 1.0) {
					float sc_x = out_x / sx;
					float sc_y = out_y / sy;
					Pixel p = GetPixelClamp((int) sc_x, (int) sc_y);
					Pixel px = GetPixelClamp((int) ceil(sc_x), (int) sc_y);
					Pixel py = GetPixelClamp((int) sc_x, (int) ceil(sc_y));
					Pixel pxy = GetPixelClamp((int) ceil(sc_x), (int) ceil(sc_y));

					float x_thingy = exp(-5*pow(ceil(sc_x) - sc_x,2));
					float y_thingy = exp(-5*pow(ceil(sc_y) - sc_y,2));

					Pixel pt = PixelLerp(p, px, x_thingy);
					Pixel pb = PixelLerp(py, pxy, x_thingy);
					out_image->GetPixel(out_x, out_y) = PixelLerp(pt, pb, y_thingy);
				}
			}
		}
	}
	return out_image;
}

Image* Image::Rotate(double angle){
	// sampling method
	return NULL;
}

void Image::Fun(){
	/* WORK HERE */
}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method){
   assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
   sampling_method = method;
}


Pixel Image::Sample (double u, double v){
   /* WORK HERE */
   return Pixel();
}