import std.math;
import std.stdio;
import std.string;

import image;
 
struct Histogram {
  int[] hist;
  int histSize;
  int below;
  int position;
  int posCount;
  int above;
}

Iplane subSample(Iplane img, int factor) {
   // already assumes size is exact factor in both x and y
   // does not handle multi-tuple images
   int row, col, index;
   int div       = factor * factor;
   int newx      = img.xsize / factor;
   int newy      = img.ysize / factor;
   double scale  = cast(double) (factor * factor);
   Iplane newImg = new Iplane("sub",newx,newy,1,img.maxval);
   for (row=0; row<img.ysize; row++) {
      for (col=0; col<img.xsize; col++) {
         index = (row / factor) * newx + (col / factor);
         newImg.img[index] += img.img[row*img.xsize+col];
      }
   }
   foreach (ref pix; newImg.img) {
      double p = floor((cast(double) pix) / scale);
      pix      = cast(int) p;
   }
   return newImg;
}

// designed to smooth background image after subsampled image expanded
Iplane blockSmooth(Iplane img, int factor) {
   int newx, newy;
   int row, col;
   int x0, x1, x2, y0, y1, y2;
   double fx0, fx1, fx2, fy0, fy1, fy2;
   double valx1y1, valx1y2, valx2y1, valx2y2;
   //writefln("BlockSmooth entered.  x=%d; y=%d; maxval=%d",img.xsize,img.ysize,img.maxval);
   newx       = img.xsize * factor;
   newy       = img.ysize * factor;
   Iplane big = new Iplane("big", newx, newy, 1, img.maxval); 
   for (row=0; row < (newy-factor); row++) {
      for (col=0; col < (newx-factor); col++) {
         // coords on old image
         x1 = col / factor;
         y1 = row / factor;
         x2 = x1 + 1;
         y2 = y1 + 1;
         valx1y1 = cast(float) img.img[y1*img.xsize+x1];
         valx1y2 = cast(float) img.img[y2*img.xsize+x1];
         valx2y1 = cast(float) img.img[y1*img.xsize+x2];
         valx2y2 = cast(float) img.img[y2*img.xsize+x2];
         // coords on new image
         fx1 = cast(double) ((col / factor) * factor);
         fy1 = cast(double) ((row / factor) * factor);
         fx2 = fx1 + cast(double) factor;
         fy2 = fy1 + cast(double) factor;
         fx0 = cast(double) col;
         fy0 = cast(double) row;
         double r1 =   ((fx2 - fx0) / (fx2 - fx1)) * valx1y1
                    +  ((fx0 - fx1) / (fx2 - fx1)) * valx2y1;
         double r2 =   ((fx2 - fx0) / (fx2 - fx1)) * valx1y2
                    +  ((fx0 - fx1) / (fx2 - fx1)) * valx2y1;
         double p  =   ((fy2 - fy0) / (fy2 - fy1)) * r1
                    +  ((fy0 - fy1) / (fy2 - fy1)) * r2;
         big.img[row*big.xsize+col] = cast(int) floor(p + 0.5);
      }
   }
   return big;
}

/*-------------------------------------------------------*/

struct MedianFilter {
   Histogram hist;
   int radius;
   int cutLevel;
   int cutVal;
   Iplane imgOut;               // medianed image
   Iplane imgIn;                // input image
   Iplane imgSub;               // background subtracted image
   private int mradius;         // used by mirror
   private int x1, x2, y1, y2;  // used by mirror, sliceout and ssMedian
  
   this(Iplane inImg, int r, int cut) {
      int bitdepth = inImg.maxval + 1;
      radius = r;
      cutLevel = cut;
      hist.histSize = bitdepth;
      hist.hist = new int[bitdepth];
      hist.below = 0;
      hist.position = 0;
      hist.posCount = 0;
      hist.above = 0;  
      imgIn      = inImg;
   } 

   Iplane doBaccor(int subsample) {
      this.ssMedian(subsample);
      return this.doSub();
   }

   Iplane doBaccor() {
      this.doMedian();
      return this.doSub();
   }

   Iplane doSub() {
      if (imgSub !is null) {
         return imgSub;
      }
      if ((imgIn is null) || (imgOut is null)) {
         throw new Exception("medianFilter.doSub() called without background");
      }

      imgSub = imgOut.blankClone();
      imgSub.img[] = imgIn.img[] - imgOut.img[];
      foreach (ref pix; imgSub.img) {
         pix = (pix >= 0)? pix : 0;
      }
      return imgSub;
   }

   Iplane doMedian() {
      if (imgOut !is null) {
         return imgOut;
      }
      int row, col, flag, index;

      flag = 0;

      this.setMirror(radius);
      imgIn   = this.mirror(imgIn);
      imgOut  = imgIn.blankClone(); 
      for (col = (radius + 1); col < (imgOut.xsize - radius - 1); col++) {
         row = radius + 1; 
         index = row * imgIn.xsize + col;
         imgOut.img[index] = initMedian(row, col);
         for (row = (radius + 2);  row < (imgIn.ysize - radius - 1); row++) {
            index = row * imgIn.xsize + col;
            imgOut.img[index] = median(row, col);
         }
      }
      imgOut = this.sliceOut(imgOut);
      imgIn  = this.sliceOut(imgIn);    // restore since needs to be intact
      return imgOut;
   }

   Iplane ssMedian(int subsample) {
      // verify subsample is possible
      // verify radius can be adjusted
      int newRadius, row, col;
      int oldX = imgIn.xsize;
      int oldY = imgIn.ysize;
      int dx   = subsample - imgIn.xsize % subsample;
      int dy   = subsample - imgIn.ysize % subsample;
      dx       = (dx==subsample)?0:dx;
      dy       = (dy==subsample)?0:dy;
      // guarantee that temp is an exact multiple of subsample
      // by padding right and bottom
      int xCorr = imgIn.xsize + dx;
      int yCorr = imgIn.ysize + dy;
      auto temp = new Iplane("",xCorr,yCorr,1,imgIn.maxval);
      // copy imgIn to temp
      for (row=0; row < imgIn.ysize; row++) {
         for (col=0; col < imgIn.xsize; col++) {
            temp.img[row*temp.xsize+col] = imgIn.img[row*imgIn.xsize+col];
         }
      }
 
      temp      = subSample(temp,subsample);
      newRadius = radius / subsample;       
      auto filt = new MedianFilter(temp,newRadius-1,cutLevel);
      temp      = filt.doMedian();
      imgOut = blockSmooth(temp,subsample);
      return imgOut;
   }

   int median (int row, int col) {
      int temp;
  
      subLine(imgIn.xsize, (col - radius), (col + radius), (row - radius - 1));
      addLine(imgIn.xsize, (col - radius), (col + radius), (row + radius));
  
      temp = walkHist();

      return temp;
   }


   int initMedian(int row, int col) {
      int x, y, i;
      int index;
      int sum, count;                      /* sum of hist stuff below median */

      for (i = 0; i < hist.histSize; i++)
        hist.hist[i] = 0;

      count = 0;
      for (y = (row - radius); y <= (row + radius); y++)
        for (x = (col - radius); x <= (col + radius); x++) {
          index = y * imgIn.xsize + x;
          hist.hist[imgIn.img[index]]++;
          count++;
        }
  
      cutVal = (count * cutLevel) / 100;
      sum = hist.position = 0;
      while (sum < cutVal) 
        sum += hist.hist[hist.position++];
  
      hist.posCount = hist.position - 1;
      hist.below = sum - hist.hist[hist.posCount];
  
      return hist.posCount;
   }


    /* 
     *  walkHist checks if hist.position still points to the median position
     *  if not, it walks the histogram in the correct direction to point to the median bin
     *  the new correct median bin index is returned
     *  This is massively faster than rescanning the histogram every time it's updated.
     */

    int walkHist()   {

        while (hist.below > cutVal) {
          hist.below -= hist.hist[--hist.posCount];
        }

        while ((hist.below + hist.hist[hist.posCount]) < cutVal) {
          hist.below += hist.hist[hist.posCount++];
        }

      /* else still okay, do nothing */
  
      return (hist.posCount);
    }


    void subLine(int xsize, int start, int end, int line) {
      int x;
      int curr;

      for (x = start; x <= end; x++) {
        curr = imgIn.img[line * xsize + x];            /* img[line,x] */
        hist.hist[curr]--;
        if (curr < hist.posCount)
          hist.below--;
      }
    }


    void addLine(int xsize, int start, int end, int line)  {
       int x;
       int curr;
  
       for (x = start; x <= end; x++) {
         curr = imgIn.img[line * xsize + x];            /* img[line,x] */
         hist.hist[curr]++;
         if (curr < hist.posCount)
           hist.below++;
       }
    }


    /* subFrom subtracts sub from modified and zeros any negative pixels
     * It does not check that the two are the same size and bit depth.
     * The programmer must guarantee that (i.e. it is for internal use not for
     * unknown files from the command line).  If external files are used, check
     * them first!!!!
     */

    Iplane subFrom(Iplane inputImg, Iplane sub) {
       int row, col, index, temp;
  
       Iplane modified = inputImg.clone();
       for (row = 0; row < modified.ysize; row++)
          for (col = 0; col < modified.xsize; col++) {
             index = row * modified.xsize + col;
             temp = modified.img[index] - sub.img[index];
             modified.img[index] = ((temp < 0) ? 0 : (temp));
          }
       return modified;
    }

   /* ---- routines to mirror sides of image ---- */

    void setMirror(int r)  {
      mradius = r;
    }

    Iplane mirror (Iplane img) {
       if (img.tuple != 1) {
          throw new Exception("Tried to mirror non-grey-scale image");
       }
       int row, col, index, index2, c, x, x2, y, y2;
       int newx, newy;
       Iplane outimg;

       newx = img.xsize + mradius * 2;
       newy = img.ysize + mradius * 2;
       this.setSliceOut(mradius,mradius,mradius+img.xsize-1,mradius+img.ysize-1);

       outimg = new Iplane("mirrored",newx, newy, 1, img.maxval);

       for (row = 0; row < img.ysize; row++)                         // copy original image
          for (col = 0; col < img.xsize; col++) {
             index = row * img.xsize + col;
             index2 = (row + mradius) * newx + (col + mradius);
             outimg.img[index2] = img.img[index];
          }

       for (row = 0; row < outimg.ysize; row++)
          for (col = 0; col < outimg.xsize; col++) {
             x2 = 0;
             y2 = 0;
             if (0 != (y2 = ytop(outimg, row, mradius))) {
	        }
             else if (0 != (y2 = ybot(outimg, row, mradius))) {
                }                
             if (0 != ( x2 = xtop(outimg, col, mradius))) {
                }
             else if (0 != (x2 = xbot(outimg, col, mradius))) {
                }


             if (x2 || y2) {
                if (x2 == 0) 
	           x2 = col;
	        if (y2 == 0)
	           y2 = row;
	        index  = row * outimg.xsize + col;
	        index2 = y2 * outimg.xsize + x2;
	        outimg.img[index] = outimg.img[index2];
             }
          }
  
      return outimg;
    }

    void setSliceOut(int xx, int yy, int xx2, int yy2) {
       x1  = xx;
       y1  = yy;
       x2 = xx2;
       y2 = yy2;
    }


    Iplane sliceOut (Iplane img) {
       int dx, dy;
       Iplane outimg;
       int row, col, indexIn, indexOut;

       dx = x2 - x1 + 1;
       dy = y2 - y1 + 1;
       outimg = new Iplane("sliced", dx, dy, 1, img.maxval);

       for (row = y1; row <= y2; row++)
          for (col = x1; col <= x2; col++) {
             indexIn = row * img.xsize + col;
             indexOut = (row - y1) * outimg.xsize + (col - x1);
             outimg.img[indexOut] = img.img[indexIn];
          }
   
       return outimg;
    }


    /*  
     *  Coordinate shift functions.
     *
     *  These functions return 0 if the coordinate doesn't need to be shifted.
     *  Otherwise, they return an offset to shift by
     *
     */

    int ytop(Iplane img, int y, int radius) {
       return ((y < radius) ? (2 * radius - y) : 0);
    }

    int ybot(Iplane img, int y, int radius) {
       return ((y>(img.ysize-radius-1)) ? (img.ysize-2*radius+(img.ysize-y)-1) : 0);
    }

    int xtop(Iplane img, int x, int radius) {
       return ((x < radius) ? (2 * radius - x) : 0); 
    }

    int xbot(Iplane img, int x, int radius) {
       return ((x>(img.xsize-radius-1)) ? (img.xsize-2*radius+(img.xsize-x)-1) : 0);
    }

}

