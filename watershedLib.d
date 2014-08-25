
/* 
 *   Vincent and Solle Watershed algorithm
 */

import std.stdio;
import std.c.stdlib;
import std.string;
import std.algorithm;
import image;

struct PointV {
   int x;
   int y;
   int val;
}

enum int Mask = -2;
enum int Wshed= 0;
enum int Init = -1;
enum int Fake = -3;

struct WQueue{
  int maxQueue;
  PointV wQueue[];
  int head = 0;
  int tail = 0;
  bool full = false;
  bool empty = true;
  
  this(int n) {
     maxQueue = n;
     wQueue   = new PointV[n];
     head = tail = 0;
     full = false; 
     empty = true;
  }

  void reset() {   // allow reuse of an old queue
     head  = tail = 0;
     full  = false;
     empty = true;
  }
  void enq(PointV pixel) {
    if (full) {
       writefln("Head = %d.  Tail = %d.", head, tail);
       throw new Exception("Watershed queue full.");
       // later, expanding the array will be an option, too lazy now
    }
    head += 1;
    head = head % maxQueue;
    wQueue[head] = pixel;
    if ((head == tail)&&(!empty)) {
      full = true;
      writefln("Head = %d.  Tail = %d.", head, tail);
      throw new Exception("Watershed queue full.");
    }
    empty = false;
    return;
  }

  PointV deq() {
     PointV temp;
     if (empty) {
        // If this happens, it is a programming error on part of the calling routine
        throw new Exception("Watershed queue empty.");
     }
     tail += 1;
     tail = tail % maxQueue;
     if (tail == head) {
         empty = true;
     } else {
         empty = false;
     }
     return wQueue[tail];
  }

  bool isEmpty() { return empty; }
} 


struct Watershedder {
   PointV[] ptArray;
   PointV[][] histArray;
   WQueue q;
   Iplane i;         // input images
   Iplane o;         // image of labeled watershed
   Iplane d;         // distances
   int min, max;
   private int currentDistance = 0;
   private int currentLabel = 0;
  
   private void hist() {
     PointV temp[];
     int first, last, n;
     int old = min = ptArray[0].val;
     max = ptArray[$-1].val;
     first = 0;
     last  = 0;
     for (n=0; n<ptArray.length; n++) {
       if (ptArray[n].val == old) {
          last = n;
       } else {
          temp = ptArray[first..(last+1)];
          histArray ~= temp;
          old        = ptArray[n].val;
          first      = n;
          last       = first;
       }
     }
     if (last >= first) {
        temp = ptArray[first..(last+1)];
        histArray ~= temp;
     }
   }
   
   this(Iplane img) {
     int row, col;
     PointV p;
     this.i = img;
     this.o = img.blankClone();
     this.o.img[] = Init;
     this.d = img.blankClone();
     this.q = WQueue(200_000);
     if (img.tuple != 1) {
        throw new Exception("Watershedded initialized with multi-or-no tuple image.");
     }
     int size = img.xsize * img.ysize;
     ptArray = new PointV[size];
     min     = this.i.maxval;
     max     = 0;
     for (row=0; row<img.ysize; row++) {
        for (col=0; col<img.xsize; col++) {
            p.val = img.img[row*img.xsize+col];
            min   = (min > p.val)? p.val : min;
            max   = (max < p.val)? p.val : max;
            p.x   = col;
            p.y   = row;
            ptArray[row*img.xsize+col] = p;
        }
     }
     sort!((a, b) { return a.val < b.val; })(ptArray);
     this.hist();
   }

   bool wORmarked(int x, int y, Iplane ip) {
      int p = ip.img[y*ip.xsize+x];
      if ((p > 0)||(p==Wshed)) {
        return true;
      } else {
        return false;
      }
   }

   // neighbor point is positive or watershed
   private bool neighbors(int x, int y, Iplane ip) {
      int xs = ip.xsize;
      int ys = ip.ysize;
      if (x==0) {
        if (wORmarked(x+1,y,ip)) { return true; }
      } else if (x==(xs-1)) {
        if (wORmarked(x-1,y,ip)) { return true; }
      } else {
        if (wORmarked(x+1,y,ip)) { return true; }
        if (wORmarked(x-1,y,ip)) { return true; }
      }
      if (y==0) {
        if (wORmarked(x,y+1,ip)) { return true; }
      } else if (y==(ys-1)) {
        if (wORmarked(x,y-1,ip)) { return true; }
      } else {
        if (wORmarked(x,y+1,ip)) { return true; }
        if (wORmarked(x,y-1,ip)) { return true; }
      }
      return false;
   }

   private PointV[] neighborVec(int x, int y, Iplane ip) {
      int n, xs, ys;
      enum Neighbors = 4;         // change if this is modified for different connectivity
      PointV[] result = new PointV[Neighbors];
      xs = ip.xsize;
      ys = ip.ysize;
      n = 0;
      if (x==0) {
        result[n] = PointV();
        result[n].x = x+1;
        result[n++].y = y;
        result.length -= 1;
      } else if (x==(xs-1)) {
        result[n] = PointV();
        result[n].x = x-1;
        result[n++].y = y;
        result.length -= 1;
      } else {
        result[n] = PointV();
        result[n].x = x+1;
        result[n++].y = y;
        result[n] = PointV();
        result[n].x = x-1;
        result[n++].y = y;
      }
      if (y==0) {
        result[n] = PointV();
        result[n].x = x;
        result[n++].y = y+1;
        result.length -= 1;
      } else if (y==(ys-1)) {
        result[n] = PointV();
        result[n].x = x;
        result[n++].y = y-1;
        result.length -= 1;
      } else {
        result[n] = PointV();
        result[n].x = x;
        result[n++].y = y+1;
        result[n] = PointV();
        result[n].x = x;
        result[n++].y = y-1;
      }
      foreach (ref PointV pp; result) {
         pp.val = ip.img[pp.y*xs+pp.x];
      }
      return result;
   }

   void floodLevel(PointV[] h) {
      PointV fake;
      fake.val = Fake;
      fake.x   = -1;
      fake.y   = -1;
      //writefln("%d pixels in h with value %d.", h.length, h[0].val);
      foreach(PointV p; h) {
         o.setPixel(p.x, p.y, Mask);
         if (neighbors(p.x, p.y, o)) {
           d.setPixel(p.x, p.y, 1);
           q.enq(p);
         } else {}
      }
      currentDistance = 1;
      q.enq(fake);
      PointV p;
      while (true) {
         p = q.deq();
         if (p.val == Fake) {
            if (q.isEmpty()) { 
               break;
            } else {
               q.enq(fake);
               currentDistance += 1;
               p = q.deq();
            }
         }
         // for all pixel p' in neighborhood p
         PointV[] pVec = neighborVec(p.x,p.y,i);
         int ip        = p.y*i.xsize+p.x;
         foreach(pp; pVec) {
            int ipp = pp.y*i.xsize+pp.x;
            // if pp belongs to an already labeled basin or watershed
            if ((d.img[ipp]<currentDistance) &&
                ((o.img[ipp]>0) || (o.img[ipp]==Wshed))) {
               if (o.img[ipp]>0) {
                  if (((o.img[ip])==Mask)||(o.img[ip]==Wshed)) {
                     o.img[ip] = o.img[ipp];
                  } else if (o.img[ip] != o.img[ipp]) {
                     o.img[ip] = Wshed;
                  } else if ((o.img[ip]==Mask)) {
                     o.img[ip] = Wshed;
                  } else {} }
                else if ((o.img[ipp]==Mask)&&(d.img[ipp]==0)) {
                     d.img[ipp] = currentDistance + 1;
                     pp.val = i.img[pp.y*i.xsize+pp.x];
                     q.enq(pp);                      
                  } else {
                    ;
                  }
               }
         }
         
      }
      foreach(PointV p1; h) {
         int ip = p1.y*i.xsize+p1.x;
         d.img[ip]  = 0;
         if (o.img[ip] == Mask) {
            currentLabel += 1;
            q.enq(p1);
            o.img[ip]     = currentLabel;
            while (!q.isEmpty()) {
               PointV p2   = q.deq();
               PointV[] p3 = neighborVec(p2.x,p2.y,i);
               foreach(pn; p3) {
                  int ipn = pn.y*i.xsize+pn.x;
                  if (o.img[ipn]==Mask) {
                     q.enq(pn);
                     o.img[ipn] = currentLabel;
                  }
               }
            } 
         }
      }
      
   }

   void process() {
      int n;
      foreach(PointV[] h; histArray) {
         if (h.length > 0) {
             //writefln("Flooding value %d.", h[0].val);
             floodLevel(h);
         }
      }
   }

   void errorCount() {
      int fake, mask, init; 
      foreach (int p; o.img) {
         switch (p) {
            case Fake: 
                fake += 1;
                break;
            case Mask:
                mask += 1;
                break;
            case Init:
                init += 1;
                break;
            default:
                ;
       }
      }
      //writeln("---- Measuring unprocessed pixels ----");
      //writefln("----\nFake = %d;  Mask = %d; Init = %d.\n---", fake, mask, init);
    }

    // write watershed lines -- needs check that img and o have same dimensions
    Iplane waterLines(Iplane img) {
       Iplane outImg = img.clone();
       int row, col, val;
       for (row=0; row<o.ysize; row++) {
          for (col=0; col<o.xsize; col++) {
             val = o.img[row*o.xsize+col];
             PointV[] nvec = this.neighborVec(col, row, o);
             foreach(PointV p; nvec) {
                if ((p.val != val)||(p.val == 0)) {
                   outImg.setPixel(p.x, p.y, 0);
                }
             }         
          }
       }
       return outImg;
    }

}
/*
void main(string args[]) {
  Iplane img;
  Watershedder w;
  if (args.length < 2) {
     writeln("watershed imageIn.pgm imageOut.pgm waterlines.pgm");
     writeln("Terminating program.");
     exit(1);
  }
     img = readPNM(args[1]);
     w   = Watershedder(img.invert());
     w.process();
     w.errorCount();
     Iplane o = w.o;
     writePNM(args[2],o);
     Iplane waterLinesImg = w.waterLines(img);
     writePNM(args[3], waterLinesImg);
}
*/