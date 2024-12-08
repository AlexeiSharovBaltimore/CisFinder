#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define gdMaxColors 256
/* Map objects */
#define OBJ_POINT    1	/* individual point */
#define OBJ_LINE     2	/* individual line */
#define OBJ_POINTF   3	/* file with multiple points */
#define OBJ_LINEF    4	/* file with multiple lines (ARC/INFO ungenerate format) */
#define OBJ_POINTFV  5	/* file with points that are colored according to their values */
#define OBJ_BITMAP   6	/* bitmap object */
#define OBJ_TEXTF    7	/* text objects in a file */
#define OBJ_FILL     8	/* fill flood the map at a particular point */
#define OBJ_SCALE    9	/* draw a space scale */
#define OBJ_TEXT    10	/* one line of text */
#define OBJ_TRANSPAR 11	/* set transparent color */
#define OBJ_POINTFV_TXT 12	/* points with values printed as text */
#define OBJ_GIFIMAGE 13	/* load a gif image */
#define OBJ_BOX      14	/* draw a box */
#define RASTER       1	/* raster input file */
#define VECTOR       2	/* drawing script as input file */
#define RIGHT  1	/* usedfor flood fill */
#define LEFT  -1	/* usedfor flood fill */
#define IDRISI	0	/* IDRISI-formated raster file */
#define ASCIIGRID	1	/* ASCIIGRID format for ARC/INFO */
#define MAXINT 65535
#define MAXCOL	5

typedef struct elem_st{	/* Element of stack for fill flood */
	int a, b, c, d;
}ELEM;

typedef struct stack_st{	/* Stack for fill flood */
	ELEM *elem;
	long num;	/* number of elements in the stack */
	long allc;	/* number of memory elements allocated */
}STACK;

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
typedef gdFont *gdFontPtr;
typedef struct param_st{
	char rastfile[80];
	int pix, n_chng;
	int drawtype;
	int backgr;
	float thresh[MAXCOL];
	int col[MAXCOL];
	int ncol;
	float xcent, ycent, dx, dy;
	float tmin;
}PARAM;

static gdImagePtr image;
static int *palette;
static int XSCREEN =  600;
static int YSCREEN =  500;
static gdFontPtr gdFontGiant;
static gdFontPtr gdFontSmall;
static gdFontPtr gdFontTiny;
static STACK stack;
static int fill_color;

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
void  check      (void *ptr);
void gdImageString(gdImagePtr im, gdFontPtr f, 
	int x, int y, unsigned char *s, int color);
void gdImageCharVertical(gdImagePtr im, gdFontPtr f, int x, int y, 
	int c, int color);
void gdImageStringVertical(gdImagePtr im, gdFontPtr f, 
	int x, int y, unsigned char *s, int color);
void   writeImage(gdImagePtr im, char *filename, PARAM *p);


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

/********************************/
void   error_message   (char *message)
/********************************/
{
	printf("ERROR: %s.\n", message);
	exit(1);
}

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
 * Define the storage for the packet accumulator
 */
//static char accum[ 256 ];

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

/*********************************************/
void  stack_push  (int a, int b, int c, int d)
/*********************************************
* Used for fill flood
*/
{
	stack.elem[stack.num].a = a;
	stack.elem[stack.num].b = b;
	stack.elem[stack.num].c = c;
	stack.elem[stack.num].d = d;
	stack.num++;
	if (stack.num >= stack.allc)
		error_message("Stack too large");
}

/*********************************************/
void  stack_pop  (int *a, int *b, int *c, int *d)
/*********************************************
* Used for fill flood
*/
{
	--stack.num;
	if (stack.num < 0)
		error_message("Empty stack\n");
	*a = stack.elem[stack.num].a;
	*b = stack.elem[stack.num].b;
	*c = stack.elem[stack.num].c;
	*d = stack.elem[stack.num].d;
}

/*********************************************/
int  getbound (int x, int y, int dir, int color)
/*********************************************
* Used for fill flood
*/
{
	int xbeg;

	if (x < 0 || x >= XSCREEN || y < 0 || y >= YSCREEN) {
		printf("%d %d\n", x, y);
		error_message ("getbound");
	}
	xbeg = x;
	while (x >= 0 && x < XSCREEN && gdImageGetPixel(image,x,y) == fill_color)
	{
		if (dir > 0 || x != xbeg)
			gdImageSetPixel(image, x, y, color);
		x += dir;
	}
	if (x != xbeg)
		x -= dir;
	return (x);
}

/*********************************************/
void tilt  (int color)
/*********************************************
* Used for fill flood
*/
{
	int updwn, left, right, row, negupdwn;
	int newl, newr, center, left1, right1;

	stack_pop (&updwn, &left, &right, &row);
	negupdwn = -updwn;
	do{
		if (row < 0 || row >= YSCREEN)
			break;
		center = left;
		while (center <= right && gdImageGetPixel(image,center,row) != fill_color)
			++center;
		if (center > right)
			break;
		newl = getbound (center, row, LEFT, color);
		newr = getbound (center, row, RIGHT, color);

		if (newl < left-1)
			stack_push(negupdwn, newl, left, row-updwn);
		else if (newl > left+1)
			stack_push(updwn, left, newl, row);
		left = newl;

		if (newr > right+1)
			stack_push(negupdwn, right, newr, row-updwn);
		else if (newr < right-1)
			stack_push(updwn, newr, right, row);
		right = newr;

		row += updwn;
	}while (newr - newl >= 0);
}

/********************************/
void  scale (float x, float y, int *ix, int *iy, PARAM *p)
/********************************/
{
	*ix = (float)XSCREEN/2 + (x - p->xcent)/p->dx;
	*iy = (float)YSCREEN/2 - (y - p->ycent)/p->dy;
}

/********************************/
void  scale_back (int ix, int iy, float *x, float *y, PARAM *p)
/********************************/
{
	*x = (ix - (float)XSCREEN/2 +0.5)*p->dx + p->xcent;
	*y = ((float)YSCREEN/2 - iy - 0.5)*p->dy + p->ycent;
}

/*********************************************/
void  fill_flood  (float x0, float y0, int color, PARAM *p)
/*********************************************
* Used for fill flood
*/
{
	int left, right, iter = 0, x, y;

	scale (x0, y0, &x, &y, p);
	if (x < 0 || x >= XSCREEN || y < 0 || y >= YSCREEN)
		return;
	fill_color = gdImageGetPixel(image,x,y);
	if (color == fill_color)
		return;
	left = getbound (x, y, LEFT, color);
	right = getbound (x, y, RIGHT, color);
	if (x >=0 && x < XSCREEN)
	{
		if (y+1 >=0 && y+1<YSCREEN)
			stack_push(1, left, right, y+1);
		if (y-1 >=0 && y-1<YSCREEN)
			stack_push(-1, left, right, y-1);
	}
	while(stack.num > 0 && iter < 1000)
	{
		tilt(color);
		++iter;
	}
}

/********************************/
void circle (int xc, int yc, int r, int color, int fill)
/********************************/
{
	int x, y, dec, iy;
	y = r;
	dec = 3-2*r;
	for (x = 0; x <= y; x++)
	{
		for (iy = -y+1; iy<y && fill >=0; ++iy)
		{
			gdImageSetPixel(image, xc+x, yc+iy, fill);
			gdImageSetPixel(image, xc+iy, yc+x, fill);
			gdImageSetPixel(image, xc-x, yc+iy, fill);
			gdImageSetPixel(image, xc+iy, yc-x, fill);
		}
		gdImageSetPixel(image, xc+x, yc+y, color);
		gdImageSetPixel(image, xc+x, yc-y, color);
		gdImageSetPixel(image, xc-x, yc+y, color);
		gdImageSetPixel(image, xc-x, yc-y, color);
		gdImageSetPixel(image, xc+y, yc+x, color);
		gdImageSetPixel(image, xc+y, yc-x, color);
		gdImageSetPixel(image, xc-y, yc+x, color);
		gdImageSetPixel(image, xc-y, yc-x, color);
		if (r >= 20) {
			gdImageSetPixel(image, xc+x, yc+y-1, color);
			gdImageSetPixel(image, xc+x, yc-y+1, color);
			gdImageSetPixel(image, xc-x, yc+y-1, color);
			gdImageSetPixel(image, xc-x, yc-y+1, color);
			gdImageSetPixel(image, xc+y-1, yc+x, color);
			gdImageSetPixel(image, xc+y-1, yc-x, color);
			gdImageSetPixel(image, xc-y+1, yc+x, color);
			gdImageSetPixel(image, xc-y+1, yc-x, color);
		}
		if ( dec >= 0 )
			dec += -4*(y--)+4;
		dec += 4*x+2;
	}
}

/********************************/
void  draw_point (float x, float y, int color, PARAM *p)
/********************************/
{
	int ix, iy, dx, dy, d0;

	if (color < 0)
		return;
	scale (x, y, &ix, &iy, p);
	ix -= p->pix/2;
	iy -= p->pix/2;
	for(dx=0;dx<p->pix;++dx)
		for(dy=0;dy<p->pix;++dy)
			gdImageSetPixel(image, ix+dx, iy+dy, palette[color]);
}

/********************************/
void  draw_circle (float x, float y, int radius, int color, int fill, PARAM *p)
/********************************/
{
	int ix, iy, dx, dy, d0;

	if (color < 0)
		return;
	scale (x, y, &ix, &iy, p);
	if (radius == 0)
		gdImageSetPixel(image, ix, iy, color);
	else
		circle(ix, iy, radius,  color, fill);
}

/********************************/
void  draw_line (float x, float y, float x1, float y1, PARAM *p,
			int color, int width, int style)
/********************************/
{
	int ix, iy, ix1, iy1;

	scale (x, y, &ix, &iy, p);
	scale (x1, y1, &ix1, &iy1, p);
	gdImageLine(image, ix, iy, ix1, iy1, color, width, style);
}

/********************************/
void  draw_box (float x, float y, float x1, float y1, PARAM *p,
			int color, int fill)
/********************************/
{
	int ix, iy, ix1, iy1, i, j;

	scale (x, y, &ix, &iy, p);
	scale (x1, y1, &ix1, &iy1, p);
	if(color >=0){
		color=palette[color];
		gdImageLine(image, ix, iy, ix, iy1, color, 1, 0);
		gdImageLine(image, ix1, iy, ix1, iy1, color, 1, 0);
		gdImageLine(image, ix, iy, ix1, iy, color, 1, 0);
		gdImageLine(image, ix, iy1, ix1, iy1, color, 1, 0);
	}
	if(fill >= 0){
		fill=palette[fill];
		if(ix1 < ix){
			int swap=ix; ix=ix1; ix1=swap;
		}
		if(iy1 < iy){
			int swap=iy; iy=iy1; iy1=swap;
		}
		for(i=ix+1; i<ix1; ++i){
			for(j=iy+1; j<iy1; ++j){
				gdImageSetPixel(image, i, j, fill);
			}
		}
	}
}


/***************************************/
void  draw_linef(int color, char *objfile, int linewid, int linestyle, PARAM *p)
/***************************************/
{
	float x, y, x1, y1;
	int i, ix, iy;
	FILE *fp;
	long num=0, done = FALSE;
	char buffer[150];

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (!done)
	{
		num = 0;
		fgets(buffer, 80, fp);
		while (fscanf(fp,"%f", &x) == 1)
		{
			fscanf(fp,"%f", &y);
			if (num)
				draw_line (x, y, x1, y1,p,palette[color], linewid, linestyle);
			x1 = x;
			y1 = y;
			++num;
		}
		if (!fgets(buffer, 80, fp))
			done = TRUE;
	}
	fclose(fp);
}

/***************************************/
void  draw_pointfv (int radius, char *objfile, PARAM *p)
/***************************************/
{
	float x, y, value;
	int i, ix, iy, color=0;
	FILE *fp;
	char buffer[100];

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fscanf(fp,"%f", &value);
		fgets(buffer, 80, fp);
		scale (x, y, &ix, &iy, p);
		if (ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
			continue;
		i = 0;
		if (color < 0)
			continue;
		if (radius >= 0)
			draw_circle(x,y,radius,palette[color], palette[color], p);
	}
	fclose(fp);
}

/***************************************/
void  draw_pointf(int color, int fill, int radius,char *objfile, PARAM *p)
/***************************************/
{
	float x, y;
	int i, ix, iy, color1, len;
	FILE *fp;
	char buffer[100];

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fgets(buffer, 80, fp);
		if (fill >=0)
			fill = palette[fill];
		if (color >= 0)
			draw_circle(x,y,radius,palette[color],fill, p);
	}
	fclose(fp);
}

/***************************************/
void  draw_textf (int color, int size, char *objfile, PARAM *p)
/***************************************/
{
	float x, y;
	int i, ix, iy, color1, len, ix1, iy1, xmin, xmax, ymin, ymax;
	FILE *fp;
	char buffer[100], *ch;
	gdFontPtr f;

	fp = fopen(objfile,"r");
	if (!fp) {
		printf("file %s not found\n", objfile);
		return;
	}
	while (fscanf(fp,"%f", &x) == 1)
	{
		fscanf(fp,"%f", &y);
		fgets(buffer, 99, fp);
		scale (x, y, &ix, &iy, p);
		i = strlen(buffer)-1;
		while (i >= 0 && (buffer[i]=='\n'||buffer[i]==' '||buffer[i]=='\t'))
			--i;
		buffer[i+1] = '\0';
		ch = buffer;
		while (ch-buffer < i && (ch[0]==' '||ch[0]=='\t'))
			++ch;
		len = strlen(ch);
		if (len == 0 || ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
			continue;
		if (size%3 == 1) f = gdFontSmall;
		else if (size%3 == 2) f = gdFontGiant;
		else f = gdFontTiny;
		if (size >= 3 && size < 6) {
			xmin = ix-len*f->w/2-2;
			xmax = ix+len*f->w/2+1;
			ymin = iy-f->h/2-1;
			ymax = iy+f->h/2+1;
			for (ix1=xmin; ix1<=xmax; ++ix1) {
				for (iy1=ymin; iy1<=ymax; ++iy1) {
					if (ix1==xmin || ix1==xmax || iy1==ymin || iy1==ymax)
						gdImageSetPixel(image, ix1, iy1, color);
					else
						gdImageSetPixel(image, ix1, iy1, 0);
				}
			}
		}
		if(size < 6){
			gdImageString(image, f, ix-len*f->w/2, iy-f->h/2, ch, color);
		}else{
			gdImageStringVertical(image, f, ix-f->h/2, iy, ch, color);
		}
	}
	fclose(fp);
}

/***************************************/
void  draw_text (int color, int size, float x, float y, char *text, PARAM *p)
/***************************************/
{
	int i, ix, iy, color1, len, ix1, iy1, xmin, xmax, ymin, ymax;
	gdFontPtr f;

	scale (x, y, &ix, &iy, p);
	len = strlen(text);
	if (len == 0 || ix < 0 || ix >= XSCREEN || iy < 0 || iy >= YSCREEN)
		return;
	if (size%3 == 1) f = gdFontSmall;
	else if (size%3 == 2) f = gdFontGiant;
	else f = gdFontTiny;
	if (size >= 3 && size < 6) {
		xmin = ix-len*f->w/2-4;
		xmax = ix+len*f->w/2+1;
		ymin = iy-f->h/2-1;
		ymax = iy+f->h/2+1;
		for (ix1=xmin; ix1<=xmax; ++ix1)
			for (iy1=ymin; iy1<=ymax; ++iy1)
				gdImageSetPixel(image, ix1, iy1, 0);
		for (ix1=xmin; ix1<=xmax; ++ix1){
			gdImageSetPixel(image, ix1, ymin, color);
			gdImageSetPixel(image, ix1, ymax, color);
		}
		for (iy1=ymin; iy1<=ymax; ++iy1){
			gdImageSetPixel(image, xmin, iy1, color);
			gdImageSetPixel(image, xmax, iy1, color);
		}
	}
	if(size < 6){
		gdImageString(image, f, ix-len*f->w/2, iy-f->h/2, text, color);
	}else{
		gdImageStringVertical(image, f, ix-f->h/2, iy, text, color);
	}
}

/***************************************/
void  set_palette  ()
/***************************************/
{
int r, g, b, i, n;
check(palette = (int*)calloc(256,sizeof(int)));
palette[0] = gdImageColorAllocate(image, 0, 0, 0);
palette[1] = gdImageColorAllocate(image, 0, 0, 127);
palette[2] = gdImageColorAllocate(image, 0, 127, 0);
palette[3] = gdImageColorAllocate(image, 0, 127, 127);
palette[4] = gdImageColorAllocate(image, 127, 0, 0);
palette[5] = gdImageColorAllocate(image, 127, 0, 127);
palette[6] = gdImageColorAllocate(image, 127, 63, 0);
palette[7] = gdImageColorAllocate(image, 127, 127, 127);
palette[8] = gdImageColorAllocate(image, 63, 63, 63);
palette[9] = gdImageColorAllocate(image, 63, 63, 255);
palette[10] = gdImageColorAllocate(image, 63, 255, 63);
palette[11] = gdImageColorAllocate(image, 63, 255, 255);
palette[12] = gdImageColorAllocate(image, 255, 63, 63);
palette[13] = gdImageColorAllocate(image, 255, 63, 255);
palette[14] = gdImageColorAllocate(image, 255, 255, 63);
palette[15] = gdImageColorAllocate(image, 255, 255, 255);
palette[16] = gdImageColorAllocate(image, 40, 160, 0);   //A
palette[17] = gdImageColorAllocate(image, 0, 0, 180); //C
palette[18] = gdImageColorAllocate(image, 220, 180, 0);   //G
palette[19] = gdImageColorAllocate(image, 170, 0, 0);   //T
/*
palette[16] = gdImageColorAllocate(image, 0, 220, 0);   //A
palette[17] = gdImageColorAllocate(image, 0, 0, 230); //C
palette[18] = gdImageColorAllocate(image, 220, 220, 0);   //G
palette[19] = gdImageColorAllocate(image, 220, 0, 0);   //T
*/
}

/***************************************/
void gdImageChar(gdImagePtr im, gdFontPtr f, int x, int y, 
	int c, int color)
/***************************************/
{
	int cx, cy;
	int px, py;
	int fline;
	cx = 0;
	cy = 0;
	if ((c < f->offset) || (c >= (f->offset + f->nchars))) {
		return;
	}
	fline = (c - f->offset) * f->h * f->w;
	for (py = y; (py < (y + f->h)); py++) {
		for (px = x; (px < (x + f->w)); px++) {
			if (f->data[fline + cy * f->w + cx]) {
				gdImageSetPixel(im, px, py, color);	
			}
			cx++;
		}
		cx = 0;
		cy++;
	}
}

/***************************************/
void gdImageString(gdImagePtr im, gdFontPtr f, 
	int x, int y, unsigned char *s, int color)
/***************************************/
{
	int i;
	int l;

	l = strlen(s);

	for (i=0; (i<l); i++) {
		gdImageChar(im, f, x, y, s[i], color);
		x += f->w;
	}
}

/***************************************/
void gdImageCharVertical(gdImagePtr im, gdFontPtr f, int x, int y, 
	int c, int color)
/***************************************/
{
	int cx, cy;
	int px, py;
	int fline;
	cx = 0;
	cy = 0;
	if ((c < f->offset) || (c >= (f->offset + f->nchars))) {
		return;
	}
	fline = (c - f->offset) * f->h * f->w;
	for (py = y-f->w; py < y; py++) {
		for (px = x; (px < (x + f->h)); px++) {
			if (f->data[fline + cx * f->w + (f->w-cy-1)]) {
				gdImageSetPixel(im, px, py, color);	
			}
			cx++;
		}
		cx = 0;
		cy++;
	}
}

/***************************************/
void gdImageStringVertical(gdImagePtr im, gdFontPtr f, 
	int x, int y, unsigned char *s, int color)
/***************************************/
{
	int i;
	int l;

	l = strlen(s);

	for (i=0; (i<l); i++) {
		gdImageCharVertical(im, f, x, y, s[i], color);
		y -= f->w;
	}
}

/*************************************/
void  modify_name  (char *nameold, char *namenew, char *modif)
/*************************************
* Replaces extention in the file name. It is used to generate
* ".doc" file names from ".img" file names, and vice versa.
*/
{
	char *ps;

	strcpy(namenew, nameold);
	ps = strrchr(namenew,'.');
	if (ps) *ps = '\0';
	strcat(namenew, modif);
}

static int
ReadColorMap(FILE *fd, int number, unsigned char (*buffer)[256])
{
       int             i;
       unsigned char   rgb[3];


       for (i = 0; i < number; ++i) {
               if (! ReadOK(fd, rgb, sizeof(rgb))) {
                       return TRUE;
               }
               buffer[CM_RED][i] = rgb[0] ;
               buffer[CM_GREEN][i] = rgb[1] ;
               buffer[CM_BLUE][i] = rgb[2] ;
       }


       return FALSE;
}

static void
ReadImage(gdImagePtr im, FILE *fd, int len, int height, unsigned char (*cmap)[256], int interlace, int ignore)
{
       unsigned char   c;      
       int             v;
       int             xpos = 0, ypos = 0, pass = 0;
       int i;
       /* Stash the color map into the image */
       for (i=0; (i<gdMaxColors); i++) {
               im->red[i] = cmap[CM_RED][i];	
               im->green[i] = cmap[CM_GREEN][i];	
               im->blue[i] = cmap[CM_BLUE][i];	
               im->open[i] = 1;
       }
       /* Many (perhaps most) of these colors will remain marked open. */
       im->colorsTotal = gdMaxColors;
       /*
       **  Initialize the Compression routines
       */
       if (! ReadOK(fd,&c,1)) {
               return; 
       }
       if (LWZReadByte(fd, TRUE, c) < 0) {
               return;
       }

       /*
       **  If this is an "uninteresting picture" ignore it.
       */
       if (ignore) {
               while (LWZReadByte(fd, FALSE, c) >= 0)
                       ;
               return;
       }

       while ((v = LWZReadByte(fd,FALSE,c)) >= 0 ) {
               /* This how we recognize which colors are actually used. */
               if (im->open[v]) {
                       im->open[v] = 0;
               }
               gdImageSetPixel(im, xpos, ypos, v);
               ++xpos;
               if (xpos == len) {
                       xpos = 0;
                       if (interlace) {
                               switch (pass) {
                               case 0:
                               case 1:
                                       ypos += 8; break;
                               case 2:
                                       ypos += 4; break;
                               case 3:
                                       ypos += 2; break;
                               }

                               if (ypos >= height) {
                                       ++pass;
                                       switch (pass) {
                                       case 1:
                                               ypos = 4; break;
                                       case 2:
                                               ypos = 2; break;
                                       case 3:
                                               ypos = 1; break;
                                       default:
                                               goto fini;
                                       }
                               }
                       } else {
                               ++ypos;
                       }
               }
               if (ypos >= height)
                       break;
       }

fini:
       if (LWZReadByte(fd,FALSE,c)>=0) {
               /* Ignore extra */
       }
}

gdImagePtr
gdImageCreateFromGif(FILE *fd)
{
       int imageNumber;
       int BitPixel;
       int ColorResolution;
       int Background;
       int AspectRatio;
       int Transparent = (-1);
       unsigned char   buf[16];
       unsigned char   c;
       unsigned char   ColorMap[3][MAXCOLORMAPSIZE];
       unsigned char   localColorMap[3][MAXCOLORMAPSIZE];
       int             imw, imh;
       int             useGlobalColormap;
       int             bitPixel;
       int             imageCount = 0;
       char            version[4];
       gdImagePtr im = 0;
       ZeroDataBlock = FALSE;

       imageNumber = 1;
       if (! ReadOK(fd,buf,6)) {
		return 0;
	}
       if (strncmp((char *)buf,"GIF",3) != 0) {
		return 0;
	}
       strncpy(version, (char *)buf + 3, 3);
       version[3] = '\0';

       if ((strcmp(version, "87a") != 0) && (strcmp(version, "89a") != 0)) {
		return 0;
	}
       if (! ReadOK(fd,buf,7)) {
		return 0;
	}
       BitPixel        = 2<<(buf[4]&0x07);
//       ColorResolution = (int) (((buf[4]&0x70)>>3)+1);
//       Background      = buf[5];
//       AspectRatio     = buf[6];

       if (BitSet(buf[4], LOCALCOLORMAP)) {    /* Global Colormap */
               if (ReadColorMap(fd, BitPixel, ColorMap)) {
			return 0;
		}
       }
       for (;;) {
               if (! ReadOK(fd,&c,1)) {
                       return 0;
               }
               if (c == ';') {         /* GIF terminator */
                       int i;
                       if (imageCount < imageNumber) {
                               return 0;
                       }
                       /* Terminator before any image was declared! */
                       if (!im) {
                              return 0;
                       }
		       /* Check for open colors at the end, so
                          we can reduce colorsTotal and ultimately
                          BitsPerPixel */
                       for (i=((im->colorsTotal-1)); (i>=0); i--) {
                               if (im->open[i]) {
                                       im->colorsTotal--;
                               } else {
                                       break;
                               }
                       } 
                       return im;
               }

               if (c == '!') {         /* Extension */
                       if (! ReadOK(fd,&c,1)) {
                               return 0;
                       }
                       DoExtension(fd, c, &Transparent);
                       continue;
               }

               if (c != ',') {         /* Not a valid start character */
                       continue;
               }

               ++imageCount;

               if (! ReadOK(fd,buf,9)) {
	               return 0;
               }

               useGlobalColormap = ! BitSet(buf[8], LOCALCOLORMAP);

               bitPixel = 1<<((buf[8]&0x07)+1);

               imw = LM_to_uint(buf[4],buf[5]);
               imh = LM_to_uint(buf[6],buf[7]);
	       if (!(im = gdImageCreate(imw, imh))) {
			 return 0;
	       }
               im->interlace = BitSet(buf[8], INTERLACE);
               if (! useGlobalColormap) {
                       if (ReadColorMap(fd, bitPixel, localColorMap)) { 
                                 return 0;
                       }
                       ReadImage(im, fd, imw, imh, localColorMap, 
                                 BitSet(buf[8], INTERLACE), 
                                 imageCount != imageNumber);
               } else {
                       ReadImage(im, fd, imw, imh,
                                 ColorMap, 
                                 BitSet(buf[8], INTERLACE), 
                                 imageCount != imageNumber);
               }
               if (Transparent != (-1)) {
                       gdImageColorTransparent(im, Transparent);
               }	   
       }
}

int gdImageColorClosest(gdImagePtr im, int r, int g, int b)
{
	int i;
	long rd, gd, bd;
	int ct = (-1);
	long mindist = 0;
	for (i=0; (i<(im->colorsTotal)); i++) {
		long dist;
		if (im->open[i]) {
			continue;
		}
		rd = (im->red[i] - r);	
		gd = (im->green[i] - g);
		bd = (im->blue[i] - b);
		dist = rd * rd + gd * gd + bd * bd;
		if ((i == 0) || (dist < mindist)) {
			mindist = dist;	
			ct = i;
		}
	}
	return ct;
}

int gdImageColorExact(gdImagePtr im, int r, int g, int b)
{
	int i;
	for (i=0; (i<(im->colorsTotal)); i++) {
		if (im->open[i]) {
			continue;
		}
		if ((im->red[i] == r) && 
			(im->green[i] == g) &&
			(im->blue[i] == b)) {
			return i;
		}
	}
	return -1;
}


/*************************************************************************/
void sortem   (long ie, float *a, int iperm, float *b,
		float *c, float *d, float *e, float *f, float *g, float *h)
/*************************************************************************/
{
	 long i, j, k, m, p, q, iring;
	 long lt[64], ut[64];
	 float ta, tb, tc, td, te, tf, tg, th, xa, xf, xg;
	 float xh, xe, xd, xc, xb;
	 int i__1;

	 /* Function Body */
	 j = ie-1;
	 m = 1;
	 i = 0;
	 iring = iperm + 1;
	 if (iperm > 7) {
		 iring = 1;
	 }

/* If this segment has more than two elements  we split it */
L10:	 if ((i__1 = j - i - 1) < 0) {
	goto L100;
	 } else if (i__1 == 0) {
	goto L90;
	 } else {
	goto L15;
	 }

/* p is the position of an arbitrary element in the segment we choose the 
* middle element. Under certain circumstances it may be advantageous 
* to choose p at random. */

L15:
	 p = (j + i) / 2;
	 ta = a[p];
	 a[p] = a[i];
	 switch (iring) {
	case 1:  goto L21;
	case 2:  goto L19;
	case 3:  goto L18;
	case 4:  goto L17;
	case 5:  goto L16;
	case 6:  goto L161;
	case 7:  goto L162;
	case 8:  goto L163;
	 }
L163:	 th = h[p];
	 h[p] = h[i];
L162:	 tg = g[p];
    g[p] = g[i];
L161:	 tf = f[p];
    f[p] = f[i];
L16:	 te = e[p];
	 e[p] = e[i];
L17:	 td = d[p];
    d[p] = d[i];
L18:	 tc = c[p];
    c[p] = c[i];
L19:	 tb = b[p];
	 b[p] = b[i];
L21: /* Start at the beginning of the segment, search for k such that a(k)>t */
    q = j;
    k = i;
L20:	 ++k;
	 if (k > q) {
	goto L60;
	 }
    if (a[k] <= ta) {
	goto L20;
    }
/* Such an element has now been found now search for a q such that a(q)<t 
* starting at the end of the segment. */
L30:  if (a[q] < ta) {
	goto L40;
	 }
    --q;
    if (q > k) {
	goto L30;
    }
    goto L50;

/* a(q) has now been found. we interchange a(q) and a(k) */

L40: xa = a[k];
    a[k] = a[q];
    a[q] = xa;
    switch (iring) {
	case 1:  goto L45;
	case 2:  goto L44;
	case 3:  goto L43;
	case 4:  goto L42;
	case 5:  goto L41;
	case 6:  goto L411;
	case 7:  goto L412;
	case 8:  goto L413;
    }
L413:     xh = h[k];
    h[k] = h[q];
	 h[q] = xh;
L412:	 xg = g[k];
    g[k] = g[q];
	 g[q] = xg;
L411:    xf = f[k];
	 f[k] = f[q];
    f[q] = xf;
L41:	 xe = e[k];
	 e[k] = e[q];
    e[q] = xe;
L42:	 xd = d[k];
    d[k] = d[q];
    d[q] = xd;
L43:    xc = c[k];
	 c[k] = c[q];
	 c[q] = xc;
L44:	 xb = b[k];
    b[k] = b[q];
	 b[q] = xb;
L45:
/* Update q and search for another pair to interchange: */
    --q;
   goto L20;
L50:    q = k - 1;
L60:
/* The upwards search has now met the downwards search: */
    a[i] = a[q];
    a[q] = ta;
   switch (iring) {
	case 1:  goto L65;
	case 2:  goto L64;
	case 3:  goto L63;
	case 4:  goto L62;
	case 5:  goto L61;
	case 6:  goto L611;
	case 7:  goto L612;
	case 8:  goto L613;
    }
L613:	 h[i] = h[q];
	 h[q] = th;
L612:    g[i] = g[q];
	 g[q] = tg;
L611:    f[i] = f[q];
    f[q] = tf;
L61:    e[i] = e[q];
	 e[q] = te;
L62:	 d[i] = d[q];
    d[q] = td;
L63:    c[i] = c[q];
    c[q] = tc;
L64:    b[i] = b[q];
    b[q] = tb;
L65:

/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
/* store the position of the largest segment in lt and ut */
    if (q << 1 <= i + j) {
	goto L70;
 }
 lt[m - 1] = i;
 ut[m - 1] = q - 1;
 i = q + 1;
 goto L80;
L70:	 lt[m - 1] = q + 1;
	 ut[m - 1] = j;
	 j = q - 1;
/* Update m and split the new smaller segment */
L80:	 ++m;
	 goto L10;

/* We arrive here if the segment has  two elements we test to see if */
/* the segment is properly ordered if not, we perform an interchange */
L90:
    if (a[i] <= a[j]) {
	goto L100;
    }
	 xa = a[i];
	 a[i] = a[j];
	 a[j] = xa;
	 switch (iring) {
	case 1:  goto L95;
	case 2:  goto L94;
	case 3:  goto L93;
	case 4:  goto L92;
	case 5:  goto L91;
	case 6:  goto L911;
	case 7:  goto L912;
	case 8:  goto L913;
    }
L913:	 xh = h[i];
	 h[i] = h[j];
	 h[j] = xh;
L912:    xg = g[i];
    g[i] = g[j];
	 g[j] = xg;
L911:	 xf = f[i];
    f[i] = f[j];
	 f[j] = xf;
L91:	 xe = e[i];
    e[i] = e[j];
    e[j] = xe;
L92:	 xd = d[i];
	 d[i] = d[j];
    d[j] = xd;
L93:	 xc = c[i];
	 c[i] = c[j];
	 c[j] = xc;
L94:    xb = b[i];
    b[i] = b[j];
	 b[j] = xb;
L95:

/* If lt and ut contain more segments to be sorted repeat process: */
L100:	 --m;
	 if (m <= 0) {
	goto L110;
	 }
	 i = lt[m - 1];
	 j = ut[m - 1];
	 goto L10;
L110:	 return;
} /* sortem_ */

/****************************/
void   check     (void *ptr)
/****************************/
{
   if (ptr == NULL)
		error_message ("Out of memory");
}

/***************************************/
int  main   (int argc, char *argv[])
/***************************************/
{
int ibp, i, j, k, n=0, motifNo=-1;
FILE *fp, *fout;
float freq[50][4];
char *buffer, motifName[100];
int color[4] = {16,17,18,19};
int *bp[4] = {NULL,NULL,NULL,NULL};

check(buffer = (char*)malloc(1500*sizeof(char)));
if(!strcmp(argv[3],"-num")){
	motifNo = atoi(argv[4]);
	if(motifNo<=0 || motifNo>10000) error_message("Wrong motif number");
}
else if (!strcmp(argv[3],"-name")) strcpy(motifName,argv[4]);
for(ibp=0; ibp<4; ++ibp){
	check(bp[ibp] = (int*)calloc(625,sizeof(int)));
	for(j=0; j<25; ++j){
		for(i=0; i<25; ++i){
			if(ibp==0){
				if(i<j/2 || i>24-j/2) continue;
				if(i>j/2+4 && i<24-j/2-4 && (j<4 || j>8)) continue;
			}else if(ibp==3){
				if(j<20 && (i<10 || i>14)) continue;
			}else if(ibp==1){
				if     (j==24 || j==0){ if(i<8 || i>16) continue; }
				else if(j==23 || j==1){ if(i<5 || i>19) continue; }
				else if(j==22 || j==2){ if (i<4 || i>21) continue; }
				else if(j==21 || j==3){ if (i<3 || i>22) continue; }
				else if(j==20 || j==4){ if (i<2 || i>22 || i>8 && i<17) continue; }
				else if(j==19 || j==5){ if (i<2 || i>23 || i>6 && i<18) continue; }
				else if(j==18 || j==6){ if (i<1 || i>23 || i>6 && i<18) continue; }
				else if(j==17 || j==7){ if (i<1 || i>6 && i<19) continue; }
				else if(j==16 || j==8){ if (i>5 && i<19) continue; }
				else if(j==15 || j==9){ if (i>5 && i<19) continue; }
				else if(j==14 || j==10){ if (i>5) continue; }
				else if(j==13 || j==11){ if (i>5) continue; }
				else if(i>5) continue;
			}else if(ibp==2){
				if(j==24 || j==0){ if(i<8 || i>16) continue; }
				else if(j==23 || j==1){ if(i<5 || i>19) continue; }
				else if(j==22 || j==2){ if (i<4 || i>21) continue; }
				else if(j==21 || j==3){ if (i<3 || i>22) continue; }
				else if(j==20 || j==4){ if (i<2 || i>22 || i>8 && i<17) continue; }
				else if(j==19 || j==5){ if (i<2 || i>23 || i>6 && i<18) continue; }
				else if(j==18 || j==6){ if (i<1 || i>23 || i>6 && i<18) continue; }
				else if(j==17 || j==7){ if (i<1 || i>6 && i<19) continue; }
				else if(j==16){ if (i>5 && i<19) continue; }
				else if(j==15){ if (i>5 && i<19) continue; }
				else if(j==14){ if (i>5) continue; }
				else if(j==13){ if (i>5) continue; }
				else if(j==11){ if (i>5 && i<15) continue; }
				else if(j==10){ if (i>5 && i<15) continue; }
				else if(j==9){ if (i>5 && i<15) continue; }
				else if(j==8){ if (i>5 && i<15) continue; }
				else if(i>5) continue;
			}
			bp[ibp][j*25+i]=1;
		}
	}
	//for(j=24; j>=0; --j){
	//	for(i=0; i<25; ++i)
	//		printf("%d",bp[ibp][j*25+i]);
	//	printf("\n");
	//}
	//printf("\n");
}
if (argc < 3)
	error_message("motiflogo inputfile outfile");
fp = fopen(argv[1],"r");
if(!fp){ printf("File %s not found\n",argv[1]); exit(0); }
int count=0;
int found = 0;
while(fgets(buffer,1499,fp)){
	//printf("%s\n",buffer);
	if(buffer[0] == '>') ++count;
	if(motifNo>=0){
		if(count==motifNo){
			found = 1;
			break;
		}
	}else{
		char *ch;
		ch = strchr(buffer,'\t');
		if(!ch) ch = strchr(buffer,'\n');
		ch[0] = '\0';
		if(!strcmp(&buffer[1],motifName)){
			found = 1;
			break;
		}
	}
}
if(!found) error_message("Motif not found!");
while(fgets(buffer,1299,fp)){
	int pos1;
	if(strlen(buffer)<3){ break; }
	if(sscanf(buffer,"%d%f%f%f%f",&pos1,&freq[n][0],&freq[n][1],&freq[n][2],&freq[n][3])<5){
		++n;
		break;
	}
	//printf("%f %f %f %f\n",freq[n][0],freq[n][1],freq[n][2],freq[n][3]);
	++n;
	if(n>=50) break;
}
XSCREEN = n*27+10;
YSCREEN = 70;

image = gdImageCreate(XSCREEN, YSCREEN);
set_palette();

/* draw background */
for(i=0; i < XSCREEN; ++i)
	for(j=0; j < YSCREEN; ++j)
		gdImageSetPixel(image, i, j, palette[15]);

float maxEntropy = 2;
float log2 = log(2.0);
for(ibp=0; ibp < n; ++ibp){
	float sum=0, info, p, index[4], data[4];
	int xstart, ystart, i1;

	info = maxEntropy;
	for(i=0; i<4; ++i) sum += freq[ibp][i];
	if(sum<=0) continue;
	for(i=0; i<4; ++i){
		freq[ibp][i] /= sum;
		index[i] = i;
		data[i] = freq[ibp][i];
		if(freq[ibp][i])
			info += freq[ibp][i]*log(freq[ibp][i])/log2;
	}
	sortem(4, data, 1, index, NULL, NULL, NULL, NULL, NULL, NULL);
	ystart = -43;
	xstart = ibp*27+5;
	for(i1=0; i1<4; ++i1){
		int ht;
		i = index[i1];
		ht = data[i1]*info*32;
		if(!ht) continue;
		for(j=0; j<25; ++j){
			for(k=0; k<ht; ++k){
				int ik, iy;
				ik = k*25/ht;
				if(!bp[i][ik*25+j] && ht>=5) continue;
				iy = 110-k;
				gdImageSetPixel(image, xstart+j, ystart+iy, color[i]);
			}
		}
		ystart -= ht;
	}
}
/* write output */
fout = fopen(argv[2], "wb");
if(!fout) error_message("Cannot write to gif file");
gdImageGif(image, fout);
fclose(fout);
gdImageDestroy(image);
return(0);
}


	


