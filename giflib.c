/*************************************
* Programmer: Alexei Sharov   11/03/2000
* Department of Entomlogy, Virginia Tech
* (540)231-7316, e-mail sharov@vt.edu
* Home page: www.ento.vt.edu/~sharov
* 
* The 'giflib.c' program has functions for generating GIF
* graphic files. It is included into program togif.c
* which is called from 'stsdecis', 'treat', and 'detailmp'
* PERL scripts
*
* Program GIFLIB.C is a part of the Slow-the-Spread
* Decision-support system. If was modified from the freeware
* software package GD.
**************************************/

#include <stdio.h>
#include <math.h>
#include "mylib.h"

#define gdMaxColors 256

typedef struct gdImageStruct {
	unsigned char ** pixels;
	int sx;
	int sy;
	int colorsTotal;
	int red[gdMaxColors];
	int green[gdMaxColors];
	int blue[gdMaxColors];
	int open[gdMaxColors];
	int transparent;
	int interlace;
} gdImage;

typedef gdImage * gdImagePtr;

typedef struct {
	/* # of characters in font */
	int nchars;
	/* First character is numbered... (usually 32 = space) */
	int offset;
	/* Character width and height */
	int w;
	int h;
	/* Font data; array of characters, one row after another.
		Easily included in code, also easily loaded from
		data files. */
	char *data;
} gdFont;

/* Text functions take these. */
typedef gdFont *gdFontPtr;

/* Functions to manipulate images. */
gdImagePtr gdImageCreate(int sx, int sy);
void gdImageDestroy(gdImagePtr im);
void gdImageSetPixel(gdImagePtr im, int x, int y, int color);
void gdImageLine(gdImagePtr im, int x1, int y1, int x2, int y2,
		int color, int width, int style);
int gdImageBoundsSafe(gdImagePtr im, int x, int y);
int gdImageColorAllocate(gdImagePtr im, int r, int g, int b);
void gdImageColorDeallocate(gdImagePtr im, int color);
void gdImageColorTransparent(gdImagePtr im, int color);
void gdImageGif(gdImagePtr im, FILE *out);
void gdImageInterlace(gdImagePtr im, int interlaceArg);

/* Macros to access information about images. READ ONLY. Changing
	these values will NOT have the desired result. */
#define gdImageSX(im) ((im)->sx)
#define gdImageSY(im) ((im)->sy)
#define gdImageColorsTotal(im) ((im)->colorsTotal)
#define gdImageRed(im, c) ((im)->red[(c)])
#define gdImageGreen(im, c) ((im)->green[(c)])
#define gdImageBlue(im, c) ((im)->blue[(c)])
#define gdImageGetTransparent(im) ((im)->transparent)
#define gdImageGetInterlaced(im) ((im)->interlace)

static unsigned int lineStyle[4] = {0xFFFF, 0xF3F3, 0x3333, 0x7272};

/***************************************/
int checkStyle (unsigned int style, unsigned int n)
/***************************************/
{
	unsigned int x=1, i;

	for(i=0;i<n%16;++i)
		x *= 2;
	if (style & x)
		return(1);
	return(0);
}

/***************************************/
int gdImageGetPixel(gdImagePtr im, int x, int y)
/***************************************/
{
	if (gdImageBoundsSafe(im, x, y)) {
		/* NOW ROW-MAJOR IN GD 1.3 */
		return im->pixels[y][x];
	} else {
		return 0;
	}
}

/***************************************/
gdImagePtr gdImageCreate(int sx, int sy)
/***************************************/
{
	int i;
	gdImagePtr im;

	check(im = (gdImage *) malloc(sizeof(gdImage)));
	check(im->pixels = (unsigned char **) malloc(sizeof(unsigned char *) * sy));
	for (i=0; (i<sy); i++) {
		check(im->pixels[i] = (unsigned char *) calloc(
			sx, sizeof(unsigned char)));
	}
	im->sx = sx;
	im->sy = sy;
	im->colorsTotal = 0;
	im->transparent = (-1);
	im->interlace = 0;
	return im;
}

/***************************************/
void gdImageDestroy(gdImagePtr im)
/***************************************/
{
	int i;
	for (i=0; (i<im->sy); i++) {
		free(im->pixels[i]);
	}	
	free(im->pixels);
	free(im);
}

/***************************************/
int gdImageColorAllocate(gdImagePtr im, int r, int g, int b)
/***************************************/
{
	int i;
	int ct = (-1);
	for (i=0; (i<(im->colorsTotal)); i++) {
		if (im->open[i]) {
			ct = i;
			break;
		}
	}	
	if (ct == (-1)) {
		ct = im->colorsTotal;
		if (ct == gdMaxColors) {
			return -1;
		}
		im->colorsTotal++;
	}
	im->red[ct] = r;
	im->green[ct] = g;
	im->blue[ct] = b;
	im->open[ct] = 0;
	return ct;
}

/***************************************/
void gdImageColorDeallocate(gdImagePtr im, int color)
/***************************************/
{
	/* Mark it open. */
	im->open[color] = 1;
}

/***************************************/
void gdImageColorTransparent(gdImagePtr im, int color)
/***************************************/
{
	im->transparent = color;
}

/***************************************/
void gdImageSetPixel(gdImagePtr im, int x, int y, int color)
/***************************************/
{
	if (gdImageBoundsSafe(im, x, y))
		im->pixels[y][x] = color;
}

/***************************************/
void gdImageSetWidePixel(gdImagePtr im, int xc, int yc, int color, int width)
/***************************************/
{
	int ix, iy, xmin, xmax, ymin, ymax;

	if (width < 2)
		gdImageSetPixel(im, xc, yc, color);
	else
	{
		xmin = xc - width/2;
		xmax = xmin + width;
		ymin = yc - width/2;
		ymax = ymin + width;
		for(ix=xmin; ix<xmax; ++ix)
			for(iy=ymin; iy<ymax; ++iy)
				gdImageSetPixel(im, ix, iy, color);
	}
}

/***************************************/
void gdImageLine(gdImagePtr im, int x1, int y1, int x2, int y2,
			int color, int width, int style)
/***************************************/
{
	int dx, dy, incr1, incr2, d, x, y, xend, yend, xdirflag, ydirflag;
	unsigned int iter=0;

	dx = abs(x2-x1);
	dy = abs(y2-y1);
	if (dy <= dx) {
		d = 2*dy - dx;
		incr1 = 2*dy;
		incr2 = 2 * (dy - dx);
		if (x1 > x2) {
			x = x2;
			y = y2;
			ydirflag = (-1);
			xend = x1;
		} else {
			x = x1;
			y = y1;
			ydirflag = 1;
			xend = x2;
		}
		if(checkStyle(lineStyle[style], iter++))
			gdImageSetWidePixel(im, x, y, color, width);
		if (((y2 - y1) * ydirflag) > 0) {
			while (x < xend) {
				x++;
				if (d <0) {
					d+=incr1;
				} else {
					y++;
					d+=incr2;
				}
				if(checkStyle(lineStyle[style], iter++))
					gdImageSetWidePixel(im, x, y, color, width);
			}
		} else {
			while (x < xend) {
				x++;
				if (d <0) {
					d+=incr1;
				} else {
					y--;
					d+=incr2;
				}
				if(checkStyle(lineStyle[style], iter++))
					gdImageSetWidePixel(im, x, y, color, width);
			}
		}
	} else {
		d = 2*dx - dy;
		incr1 = 2*dx;
		incr2 = 2 * (dx - dy);
		if (y1 > y2) {
			y = y2;
			x = x2;
			yend = y1;
			xdirflag = (-1);
		} else {
			y = y1;
			x = x1;
			yend = y2;
			xdirflag = 1;
		}
		if(checkStyle(lineStyle[style], iter++))
			gdImageSetWidePixel(im, x, y, color, width);
		if (((x2 - x1) * xdirflag) > 0) {
			while (y < yend) {
				y++;
				if (d <0) {
					d+=incr1;
				} else {
					x++;
					d+=incr2;
				}
				if(checkStyle(lineStyle[style], iter++))
					gdImageSetWidePixel(im, x, y, color, width);
			}
		} else {
			while (y < yend) {
				y++;
				if (d <0) {
					d+=incr1;
				} else {
					x--;
					d+=incr2;
				}
				if(checkStyle(lineStyle[style], iter++))
					gdImageSetWidePixel(im, x, y, color, width);
			}
		}
	}
}

/***************************************/
int gdImageBoundsSafe(gdImagePtr im, int x, int y)
/***************************************/
{
	return (!(((y < 0) || (y >= im->sy)) ||
		((x < 0) || (x >= im->sx))));
}

/* Code drawn from ppmtogif.c, from the pbmplus package
**
** Based on GIFENCOD by David Rowley <mgardi@watdscu.waterloo.edu>. A
** Lempel-Zim compression based on "compress".
**
** Modified by Marcel Wijkstra <wijkstra@fwi.uva.nl>
**
** Copyright (C) 1989 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
**
** The Graphics Interchange Format(c) is the Copyright property of
** CompuServe Incorporated.  GIF(sm) is a Service Mark property of
** CompuServe Incorporated.
*
*  Heavily modified by Mouse, 1998-02-12.  
*  Remove LZW compression.
*  Added miGIF run length compression.
*
*/

/*
 * a code_int must be able to hold 2**GIFBITS values of type int, and also -1
 */
typedef int code_int;

static int colorstobpp(int colors);
static void BumpPixel (void);
static int GIFNextPixel (gdImagePtr im);
static void GIFEncode (FILE *fp, int GWidth, int GHeight, int GInterlace, int Background, int Transparent, int BitsPerPixel, int *Red, int *Green, int *Blue, gdImagePtr im);
static void Putword (int w, FILE *fp);
static void compress (int, FILE *, gdImagePtr, int);
static void output (code_int code);
static void char_init (void);
static void char_out (int c);
/* Allows for reuse */
static void init_statics(void);

void gdImageGif(gdImagePtr im, FILE *out)
{
	int interlace, transparent, BitsPerPixel;
	interlace = im->interlace;
	transparent = im->transparent;

	BitsPerPixel = colorstobpp(im->colorsTotal);
	/* Clear any old values in statics strewn through the GIF code */
	init_statics();
	/* All set, let's do it. */
	GIFEncode(
		out, im->sx, im->sy, interlace, 0, transparent, BitsPerPixel,
		im->red, im->green, im->blue, im);
}

static int
colorstobpp(int colors)
{
    int bpp = 0;

    if ( colors <= 2 )
        bpp = 1;
    else if ( colors <= 4 )
        bpp = 2;
    else if ( colors <= 8 )
		  bpp = 3;
    else if ( colors <= 16 )
        bpp = 4;
    else if ( colors <= 32 )
        bpp = 5;
    else if ( colors <= 64 )
        bpp = 6;
    else if ( colors <= 128 )
        bpp = 7;
    else if ( colors <= 256 )
        bpp = 8;
    return bpp;
    }

/*****************************************************************************
 *
 * GIFENCODE.C    - GIF Image compression interface
 *
 * GIFEncode( FName, GHeight, GWidth, GInterlace, Background, Transparent,
 *            BitsPerPixel, Red, Green, Blue, gdImagePtr )
 *
 *****************************************************************************/

#define TRUE 1
#define FALSE 0

static int Width, Height;
static int curx, cury;
static long CountDown;
static int Pass = 0;
static int Interlace;

/*
 * Bump the 'curx' and 'cury' to point to the next pixel
 */
static void
BumpPixel(void)
{
        /*
			* Bump the current X position
         */
        ++curx;

        /*
         * If we are at the end of a scan line, set curx back to the beginning
         * If we are interlaced, bump the cury to the appropriate spot,
         * otherwise, just increment it.
         */
        if( curx == Width ) {
                curx = 0;

                if( !Interlace )
								++cury;
                else {
                     switch( Pass ) {

                       case 0:
                          cury += 8;
                          if( cury >= Height ) {
                                ++Pass;
                                cury = 4;
                          }
                          break;

                       case 1:
								  cury += 8;
                          if( cury >= Height ) {
                                ++Pass;
                                cury = 2;
                          }
                          break;

                       case 2:
                          cury += 4;
                          if( cury >= Height ) {
                             ++Pass;
                             cury = 1;
                          }
								  break;

                       case 3:
                          cury += 2;
                          break;
                        }
                }
        }
}

/*
 * Return the next pixel from the image
 */
static int
GIFNextPixel(gdImagePtr im)
{
        int r;

        if( CountDown == 0 )
                return EOF;

        --CountDown;

        r = gdImageGetPixel(im, curx, cury);

        BumpPixel();

        return r;
}

/* public */

static void
GIFEncode(FILE *fp, int GWidth, int GHeight, int GInterlace, int Background, int Transparent, int BitsPerPixel, int *Red, int *Green, int *Blue, gdImagePtr im)
{
        int B;
        int RWidth, RHeight;
        int LeftOfs, TopOfs;
        int Resolution;
		  int ColorMapSize;
        int InitCodeSize;
        int i;

        Interlace = GInterlace;

        ColorMapSize = 1 << BitsPerPixel;

        RWidth = Width = GWidth;
        RHeight = Height = GHeight;
        LeftOfs = TopOfs = 0;

        Resolution = BitsPerPixel;

        /*
         * Calculate number of bits we are expecting
         */
        CountDown = (long)Width * (long)Height;

        /*
         * Indicate which pass we are on (if interlace)
         */
        Pass = 0;

        /*
         * The initial code size
			*/
        if( BitsPerPixel <= 1 )
                InitCodeSize = 2;
        else
                InitCodeSize = BitsPerPixel;

        /*
         * Set up the current x and y position
         */
        curx = cury = 0;

        /*
         * Write the Magic header
			*/
        fwrite( Transparent < 0 ? "GIF87a" : "GIF89a", 1, 6, fp );

        /*
         * Write out the screen width and height
         */
        Putword( RWidth, fp );
        Putword( RHeight, fp );

        /*
         * Indicate that there is a global colour map
         */
        B = 0x80;       /* Yes, there is a color map */

        /*
         * OR in the resolution
         */
        B |= (Resolution - 1) << 4;

        /*
         * OR in the Bits per Pixel
         */
        B |= (BitsPerPixel - 1);

        /*
         * Write it out
			*/
        fputc( B, fp );

        /*
         * Write out the Background colour
         */
        fputc( Background, fp );

        /*
         * Byte of 0's (future expansion)
         */
        fputc( 0, fp );

		  /*
         * Write out the Global Colour Map
         */
        for( i=0; i<ColorMapSize; ++i ) {
                fputc( Red[i], fp );
                fputc( Green[i], fp );
                fputc( Blue[i], fp );
        }

	/*
	 * Write out extension for transparent colour index, if necessary.
	 */
	if ( Transparent >= 0 ) {
		 fputc( '!', fp );
		 fputc( 0xf9, fp );
		 fputc( 4, fp );
		 fputc( 1, fp );
		 fputc( 0, fp );
		 fputc( 0, fp );
		 fputc( (unsigned char) Transparent, fp );
		 fputc( 0, fp );
	}

		  /*
			* Write an Image separator
			*/
		  fputc( ',', fp );

		  /*
			* Write the Image header
			*/

		  Putword( LeftOfs, fp );
		  Putword( TopOfs, fp );
		  Putword( Width, fp );
		  Putword( Height, fp );

		  /*
			* Write out whether or not the image is interlaced
			*/
		  if( Interlace )
					 fputc( 0x40, fp );
		  else
					 fputc( 0x00, fp );

		  /*
			* Write out the initial code size
			*/
		  fputc( InitCodeSize, fp );

		  /*
			* Go and actually compress the data
			*/
		  compress( InitCodeSize+1, fp, im, Background );

		  /*
			* Write out a Zero-length packet (to end the series)
			*/
		  fputc( 0, fp );

		  /*
			* Write the GIF file terminator
			*/
		  fputc( ';', fp );
}

/*
 * Write out a word to the GIF file
 */
static void
Putword(int w, FILE *fp)
{
		  fputc( w & 0xff, fp );
		  fputc( (w / 256) & 0xff, fp );
}

#define GIFBITS 12

/*-----------------------------------------------------------------------
 *
 * miGIF Compression - mouse and ivo's GIF-compatible compression
 *
 *          -run length encoding compression routines-
 *
 * Copyright (C) 1998 Hutchison Avenue Software Corporation
 *               http://www.hasc.com
 *               info@hasc.com
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose and without fee is hereby granted, provided
 * that the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  This software is provided "AS IS." The Hutchison Avenue 
 * Software Corporation disclaims all warranties, either express or implied, 
 * including but not limited to implied warranties of merchantability and 
 * fitness for a particular purpose, with respect to this code and accompanying
 * documentation. 
 * 
 * The miGIF compression routines do not, strictly speaking, generate files 
 * conforming to the GIF spec, since the image data is not LZW-compressed 
 * (this is the point: in order to avoid transgression of the Unisys patent 
 * on the LZW algorithm.)  However, miGIF generates data streams that any
 * reasonably sane LZW decompresser will decompress to what we want.
 *
 * miGIF compression uses run length encoding. It compresses horizontal runs 
 * of pixels of the same color. This type of compression gives good results
 * on images with many runs, for example images with lines, text and solid 
 * shapes on a solid-colored background. It gives little or no compression 
 * on images with few runs, for example digital or scanned photos.
 *
 *                               der Mouse
 *                      mouse@rodents.montreal.qc.ca
 *            7D C8 61 52 5D E7 2D 39  4E F1 31 3E E8 B3 27 4B
 *
 *                             ivo@hasc.com
 *
 * The Graphics Interchange Format(c) is the Copyright property of
 * CompuServe Incorporated.  GIF(sm) is a Service Mark property of
 * CompuServe Incorporated.
 *
 */

static int rl_pixel;
static int rl_basecode;
static int rl_count;
static int rl_table_pixel;
static int rl_table_max;
static int just_cleared;
static int out_bits;
static int out_bits_init;
static int out_count;
static int out_bump;
static int out_bump_init;
static int out_clear;
static int out_clear_init;
static int max_ocodes;
static int code_clear;
static int code_eof;
static unsigned int obuf;
static int obits;
static FILE *ofile;
static unsigned char oblock[256];
static int oblen;

/* Used only when debugging GIF compression code */
/* #define DEBUGGING_ENVARS */

#ifdef DEBUGGING_ENVARS

static int verbose_set = 0;
static int verbose;
#define VERBOSE (verbose_set?verbose:set_verbose())

static int set_verbose(void)
{
 verbose = !!getenv("GIF_VERBOSE");
 verbose_set = 1;
 return(verbose);
}

#else

#define VERBOSE 0

#endif


static const char *binformat(unsigned int v, int nbits)
{
 static char bufs[8][64];
 static int bhand = 0;
 unsigned int bit;
 int bno;
 char *bp;

 bhand --;
 if (bhand < 0) bhand = (sizeof(bufs)/sizeof(bufs[0]))-1;
 bp = &bufs[bhand][0];
 for (bno=nbits-1,bit=1U<<bno;bno>=0;bno--,bit>>=1)
  { *bp++ = (v & bit) ? '1' : '0';
    if (((bno&3) == 0) && (bno != 0)) *bp++ = '.';
  }
 *bp = '\0';
 return(&bufs[bhand][0]);
}

static void write_block(void)
{
 int i;

 fputc(oblen,ofile);
 fwrite(&oblock[0],1,oblen,ofile);
 oblen = 0;
}

static void block_out(unsigned char c)
{
 oblock[oblen++] = c;
 if (oblen >= 255) write_block();
}

static void block_flush(void)
{
 if (oblen > 0) write_block();
}

static void output(int val)
{
 obuf |= val << obits;
 obits += out_bits;
 while (obits >= 8)
  { block_out(obuf&0xff);
    obuf >>= 8;
    obits -= 8;
  }
}

static void output_flush(void)
{
 if (obits > 0) block_out(obuf);
 block_flush();
}

static void did_clear(void)
{
 out_bits = out_bits_init;
 out_bump = out_bump_init;
 out_clear = out_clear_init;
 out_count = 0;
 rl_table_max = 0;
 just_cleared = 1;
}

static void output_plain(int c)
{
 just_cleared = 0;
 output(c);
 out_count ++;
 if (out_count >= out_bump)
  { out_bits ++;
    out_bump += 1 << (out_bits - 1);
  }
 if (out_count >= out_clear)
  { output(code_clear);
    did_clear();
  }
}

/*static unsigned int isqrt(unsigned int) __attribute__((__const__)); */
static unsigned int isqrt(unsigned int x)
{
 unsigned int r;
 unsigned int v;

 if (x < 2) return(x);
 for (v=x,r=1;v;v>>=2,r<<=1) ;
 while (1)
  { v = ((x / r) + r) / 2;
    if ((v == r) || (v == r+1)) return(r);
    r = v;
  }
  return(r);
}

static unsigned int compute_triangle_count(unsigned int count, unsigned int nrepcodes)
{
 unsigned int perrep;
 unsigned int cost;

 cost = 0;
 perrep = (nrepcodes * (nrepcodes+1)) / 2;
 while (count >= perrep)
  { cost += nrepcodes;
    count -= perrep;
  }
 if (count > 0)
  { unsigned int n;
    n = isqrt(count);
    while ((n*(n+1)) >= 2*count) n --;
    while ((n*(n+1)) < 2*count) n ++;
    cost += n;
  }
 return(cost);
}

static void max_out_clear(void)
{
 out_clear = max_ocodes;
}

static void reset_out_clear(void)
{
 out_clear = out_clear_init;
 if (out_count >= out_clear)
  { output(code_clear);
    did_clear();
  }
}

static void rl_flush_fromclear(int count)
{
 int n;

 max_out_clear();
 rl_table_pixel = rl_pixel;
 n = 1;
 while (count > 0)
  { if (n == 1)
     { rl_table_max = 1;
       output_plain(rl_pixel);
       count --;
     }
    else if (count >= n)
     { rl_table_max = n;
       output_plain(rl_basecode+n-2);
       count -= n;
     }
    else if (count == 1)
     { rl_table_max ++;
       output_plain(rl_pixel);
       count = 0;
     }
    else
     { rl_table_max ++;
       output_plain(rl_basecode+count-2);
       count = 0;
     }
    if (out_count == 0) n = 1; else n ++;
  }
 reset_out_clear();
}

static void rl_flush_clearorrep(int count)
{
 int withclr;

 withclr = 1 + compute_triangle_count(count,max_ocodes);
 if (withclr < count)
  { output(code_clear);
    did_clear();
    rl_flush_fromclear(count);
  }
 else
  { for (;count>0;count--) output_plain(rl_pixel);
  }
}

static void rl_flush_withtable(int count)
{
 int repmax;
 int repleft;
 int leftover;

 repmax = count / rl_table_max;
 leftover = count % rl_table_max;
 repleft = (leftover ? 1 : 0);
 if (out_count+repmax+repleft > max_ocodes)
  { repmax = max_ocodes - out_count;
    leftover = count - (repmax * rl_table_max);
    repleft = 1 + compute_triangle_count(leftover,max_ocodes);
  }
 if (1+compute_triangle_count(count,max_ocodes) < repmax+repleft)
  { output(code_clear);
    did_clear();
    rl_flush_fromclear(count);
    return;
  }
 max_out_clear();
 for (;repmax>0;repmax--) output_plain(rl_basecode+rl_table_max-2);
 if (leftover)
  { if (just_cleared)
     { rl_flush_fromclear(leftover);
     }
    else if (leftover == 1)
     { output_plain(rl_pixel);
     }
    else
     { output_plain(rl_basecode+leftover-2);
     }
  }
 reset_out_clear();
}

static void rl_flush(void)
{
 int table_reps;
 int table_extra;

 if (rl_count == 1)
  { output_plain(rl_pixel);
    rl_count = 0;
    return;
  }
 if (just_cleared)
  { rl_flush_fromclear(rl_count);
  }
 else if ((rl_table_max < 2) || (rl_table_pixel != rl_pixel))
  { rl_flush_clearorrep(rl_count);
  }
 else
  { rl_flush_withtable(rl_count);
  }
 rl_count = 0;
}

static void compress(int init_bits, FILE *outfile, gdImagePtr im, int background)
{
 int c;

 ofile = outfile;
 obuf = 0;
 obits = 0;
 oblen = 0;
 code_clear = 1 << (init_bits - 1);
 code_eof = code_clear + 1;
 rl_basecode = code_eof + 1;
 out_bump_init = (1 << (init_bits - 1)) - 1;
 /* for images with a lot of runs, making out_clear_init larger will
    give better compression. */ 
 out_clear_init = (init_bits <= 3) ? 9 : (out_bump_init-1);
#ifdef DEBUGGING_ENVARS
  { const char *ocienv;
    ocienv = getenv("GIF_OUT_CLEAR_INIT");
    if (ocienv)
     { out_clear_init = atoi(ocienv);
       if (VERBOSE) printf("[overriding out_clear_init to %d]\n",out_clear_init);
     }
  }
#endif
 out_bits_init = init_bits;
 max_ocodes = (1 << GIFBITS) - ((1 << (out_bits_init - 1)) + 3);
 did_clear();
 output(code_clear);
 rl_count = 0;
 while (1)
  { c = GIFNextPixel(im);
	 if ((rl_count > 0) && (c != rl_pixel || rl_count==MAXINT))
		rl_flush();
	 if (c == EOF) break;
	 if (rl_pixel == c)
	  { rl_count ++;
	  }
	 else
	  { rl_pixel = c;
		 rl_count = 1;
	  }
  }
 output(code_eof);
 output_flush();
}

/*-----------------------------------------------------------------------
 *
 * End of miGIF section  - See copyright notice at start of section.
 *
/*-----------------------------------------------------------------------


/******************************************************************************
 *
 * GIF Specific routines
 *
 ******************************************************************************/

/*
 * Number of characters so far in this 'packet'
 */
static int a_count;

/*
 * Set up the 'byte output' routine
 */
static void
char_init(void)
{
        a_count = 0;
}

/*
 * Define the storage for the packet accumulator
 */
static char accum[ 256 ];

static void init_statics(void) {
	/* Some of these are properly initialized later. What I'm doing
		here is making sure code that depends on C's initialization
		of statics doesn't break when the code gets called more
		than once. */
	Width = 0;
	Height = 0;
	curx = 0;
	cury = 0;
	CountDown = 0;
	Pass = 0;
	Interlace = 0;
	a_count = 0;
}


/* +-------------------------------------------------------------------+ */
/* | Copyright 1990, 1991, 1993, David Koblas.  (koblas@netcom.com)    | */
/* |   Permission to use, copy, modify, and distribute this software   | */
/* |   and its documentation for any purpose and without fee is hereby | */
/* |   granted, provided that the above copyright notice appear in all | */
/* |   copies and that both that copyright notice and this permission  | */
/* |   notice appear in supporting documentation.  This software is    | */
/* |   provided "as is" without express or implied warranty.           | */
/* +-------------------------------------------------------------------+ */


#define        MAXCOLORMAPSIZE         256

#define        TRUE    1
#define        FALSE   0

#define CM_RED         0
#define CM_GREEN       1
#define CM_BLUE                2

#define        MAX_LWZ_BITS            12

#define INTERLACE              0x40
#define LOCALCOLORMAP  0x80
#define BitSet(byte, bit)      (((byte) & (bit)) == (bit))

#define        ReadOK(file,buffer,len) (fread(buffer, len, 1, file) != 0)

//#define LM_to_uint(a,b)                        (((b)<<8)|(a))

/* We may eventually want to use this information, but def it out for now */
#if 0
static struct {
       unsigned int    Width;
       unsigned int    Height;
       unsigned char   ColorMap[3][MAXCOLORMAPSIZE];
       unsigned int    BitPixel;
       unsigned int    ColorResolution;
       unsigned int    Background;
       unsigned int    AspectRatio;
} GifScreen;
#endif

static struct {
       int     transparent;
       int     delayTime;
       int     inputFlag;
       int     disposal;
} Gif89 = { -1, -1, -1, 0 };

static int ReadColorMap (FILE *fd, int number, unsigned char (*buffer)[256]);
static int DoExtension (FILE *fd, int label, int *Transparent);
static int GetDataBlock (FILE *fd, unsigned char *buf);
static int GetCode (FILE *fd, int code_size, int flag);
static int LWZReadByte (FILE *fd, int flag, int input_code_size);
static void ReadImage (gdImagePtr im, FILE *fd, int len, int height, unsigned char (*cmap)[256], int interlace, int ignore);
int ZeroDataBlock;

unsigned int LM_to_uint(unsigned char a, unsigned char b){
	unsigned int i;
	i = b*256+a;
	return(i);
}

static int
DoExtension(FILE *fd, int label, int *Transparent)
{
		 static unsigned char     buf[256];

		 switch (label) {
		 case 0xf9:              /* Graphic Control Extension */
					(void) GetDataBlock(fd, (unsigned char*) buf);
					Gif89.disposal    = (buf[0] >> 2) & 0x7;
					Gif89.inputFlag   = (buf[0] >> 1) & 0x1;
					Gif89.delayTime   = LM_to_uint(buf[1],buf[2]);
					if ((buf[0] & 0x1) != 0)
							  *Transparent = buf[3];

					while (GetDataBlock(fd, (unsigned char*) buf) != 0)
							  ;
					return FALSE;
		 default:
					break;
		 }
		 while (GetDataBlock(fd, (unsigned char*) buf) != 0)
               ;

       return FALSE;
}

static int
GetDataBlock_(FILE *fd, unsigned char *buf)
{
       unsigned char   count;

       if (! ReadOK(fd,&count,1)) {
               return -1;
       }

       ZeroDataBlock = count == 0;

       if ((count != 0) && (! ReadOK(fd, buf, count))) {
               return -1;
       }

       return count;
}

static int
GetDataBlock(FILE *fd, unsigned char *buf)
{
 int rv;
 int i;

 rv = GetDataBlock_(fd,buf);
 return(rv);
}

static int
GetCode_(FILE *fd, int code_size, int flag)
{
       static unsigned char    buf[280];
       static int              curbit, lastbit, done, last_byte;
       int                     i, j, ret;
       unsigned char           count;

       if (flag) {
               curbit = 0;
               lastbit = 0;
               done = FALSE;
               return 0;
       }

       if ( (curbit+code_size) >= lastbit) {
               if (done) {
                       if (curbit >= lastbit) {
                                /* Oh well */
                       }                        
                       return -1;
               }
               buf[0] = buf[last_byte-2];
               buf[1] = buf[last_byte-1];

               if ((count = GetDataBlock(fd, &buf[2])) == 0)
                       done = TRUE;

               last_byte = 2 + count;
               curbit = (curbit - lastbit) + 16;
               lastbit = (2+count)*8 ;
       }

       ret = 0;
       for (i = curbit, j = 0; j < code_size; ++i, ++j)
               ret |= ((buf[ i / 8 ] & (1 << (i % 8))) != 0) << j;

       curbit += code_size;
       return ret;
}

static int
GetCode(FILE *fd, int code_size, int flag)
{
 int rv;

 rv = GetCode_(fd,code_size,flag);
 return(rv);
}

static int
LWZReadByte_(FILE *fd, int flag, int input_code_size)
{
       static int      fresh = FALSE;
       int             code, incode;
       static int      code_size, set_code_size;
       static int      max_code, max_code_size;
       static int      firstcode, oldcode;
       static int      clear_code, end_code;
       static int      table[2][(1<< MAX_LWZ_BITS)];
       static int      stack[(1<<(MAX_LWZ_BITS))*2], *sp;
       register int    i;

       if (flag) {
               set_code_size = input_code_size;
               code_size = set_code_size+1;
               clear_code = 1 << set_code_size ;
               end_code = clear_code + 1;
               max_code_size = 2*clear_code;
               max_code = clear_code+2;

               GetCode(fd, 0, TRUE);
               
               fresh = TRUE;

               for (i = 0; i < clear_code; ++i) {
                       table[0][i] = 0;
                       table[1][i] = i;
               }
               for (; i < (1<<MAX_LWZ_BITS); ++i)
                       table[0][i] = table[1][0] = 0;

               sp = stack;

               return 0;
       } else if (fresh) {
               fresh = FALSE;
               do {
                       firstcode = oldcode =
                               GetCode(fd, code_size, FALSE);
               } while (firstcode == clear_code);
               return firstcode;
       }

       if (sp > stack)
               return *--sp;

       while ((code = GetCode(fd, code_size, FALSE)) >= 0) {
               if (code == clear_code) {
                       for (i = 0; i < clear_code; ++i) {
                               table[0][i] = 0;
                               table[1][i] = i;
                       }
                       for (; i < (1<<MAX_LWZ_BITS); ++i)
                               table[0][i] = table[1][i] = 0;
                       code_size = set_code_size+1;
                       max_code_size = 2*clear_code;
                       max_code = clear_code+2;
                       sp = stack;
                       firstcode = oldcode =
                                       GetCode(fd, code_size, FALSE);
                       return firstcode;
               } else if (code == end_code) {
                       int             count;
                       unsigned char   buf[260];

                       if (ZeroDataBlock)
                               return -2;

                       while ((count = GetDataBlock(fd, buf)) > 0)
                               ;

                       if (count != 0)
                       return -2;
               }

               incode = code;

               if (code >= max_code) {
                       *sp++ = firstcode;
                       code = oldcode;
               }

               while (code >= clear_code) {
                       *sp++ = table[1][code];
                       if (code == table[0][code]) {
                               /* Oh well */
                       }
                       code = table[0][code];
               }

               *sp++ = firstcode = table[1][code];

               if ((code = max_code) <(1<<MAX_LWZ_BITS)) {
                       table[0][code] = oldcode;
                       table[1][code] = firstcode;
                       ++max_code;
                       if ((max_code >= max_code_size) &&
                               (max_code_size < (1<<MAX_LWZ_BITS))) {
                               max_code_size *= 2;
                               ++code_size;
                       }
               }

               oldcode = incode;

               if (sp > stack)
                       return *--sp;
       }
       return code;
}

static int
LWZReadByte(FILE *fd, int flag, int input_code_size)
{
 int rv;

 rv = LWZReadByte_(fd,flag,input_code_size);
 return(rv);
}

	
