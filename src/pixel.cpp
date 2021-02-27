#include "pixel.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>



/**
 * Component
 **/
Component ComponentRandom(void)
{
    return rand() % 256;
}


Component ComponentScale(Component c, double f)
{
    return ComponentClamp((int) floor(c * f + 0.5));
}


Component ComponentLerp(Component c, Component d, double t)
{
    return ComponentClamp((int) floor((1.0 - t) * c + t * d + 0.5));
}



/// Pixel
void Pixel::SetClamp (double r_, double g_, double b_)
{
    r = ComponentClamp((int)r_);
    g = ComponentClamp((int)g_);
    b = ComponentClamp((int)b_);
}

void Pixel::SetClamp (double r_, double g_, double b_, double a_)
{
    r = ComponentClamp((int)r_);
    g = ComponentClamp((int)g_);
    b = ComponentClamp((int)b_);
    a = ComponentClamp((int)a_);
}

Pixel PixelRandom(void)
{
    return Pixel(
        ComponentRandom(),
        ComponentRandom(),
        ComponentRandom(),
        ComponentRandom());
}


Pixel operator+ (const Pixel& p, const Pixel& q)
{
    return Pixel(
        ComponentClamp(p.r + q.r),
        ComponentClamp(p.g + q.g),
        ComponentClamp(p.b + q.b),
        ComponentClamp(p.a + q.a));
}


Pixel operator* (const Pixel& p, const Pixel& q)
{
    return Pixel(
        ComponentClamp(p.r * q.r),
        ComponentClamp(p.g * q.g),
        ComponentClamp(p.b * q.b),
        ComponentClamp(p.a * q.a));
}


Pixel operator* (const Pixel& p, double f)
{
    return Pixel(
        ComponentScale(p.r, f),
        ComponentScale(p.g, f),
        ComponentScale(p.b, f),
        ComponentScale(p.a, f));
}


Pixel PixelLerp (const Pixel& p, const Pixel& q, double t)
{
    return Pixel(
        ComponentLerp(p.r, q.r, t),
        ComponentLerp(p.g, q.g, t),
        ComponentLerp(p.b, q.b, t),
        ComponentLerp(p.a, q.a, t));
}

Pixel PixelQuant( const Pixel &p, int nbits)
{
	int shift = 8-nbits;
	float mult = 255.0/(255 >> shift);
	int new_r, new_g, new_b;
	new_r = (p.r >> shift);
	new_g = (p.g >> shift);
	new_b = (p.b >> shift);

	Pixel ret;
	ret.SetClamp(new_r*mult , new_g*mult , new_b*mult );
	return ret;
}

Pixel Grayscale(float lum) {
	return Pixel(lum, lum, lum);
}

float Pixel::Hue() {
	float fr = r / 255.0;
	float fg = g / 255.0;
	float fb = b / 255.0;

	float min = fmin(fr, fmin(fg, fb));
	float max = fmax(fr, fmax(fg, fb));
	
	if (min == max) { return 0; }
	
	float hue;

	if (fr == max) {
		hue = 0.0 + (fg - fb)/(max - min);
	} else if (fg == max) {
		hue = 2.0 + (fb - fr)/(max - min);
	} else {
		hue = 4.0 + (fr - fg)/(max - min);
	}
	hue /= 6.0;
	return hue;
}

float Pixel::Saturation() {
	float fr = r / 255.0;
	float fg = g / 255.0;
	float fb = b / 255.0;

	float min = fmin(fr, fmin(fg, fb));
	float max = fmax(fr, fmax(fg, fb));
	if (Luminance() <= 0.5) {
		return (max-min)/(max+min);
	} else {
		return (max-min)/(2.0-max-min);
	}
}

float Pixel::Luminance () {
    return ((r * 76 + g * 150 + b * 29) >> 8)/255.0;
}