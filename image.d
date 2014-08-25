import std.stdio;
import std.file;
import std.string;
import std.math;
import std.conv;
import std.ascii; 
            
enum ReadState {
  Begin,
  Whitespace,
  Digit,
  Newline,
  Comment,
  Other, 
  Finished
}


class Iplane {
   string name;
   int    xsize;
   int    ysize;
   int    tuple;
   int    maxval;
   private int _padx=0;
   private int _pady=0;
   private int    _size;
   int[]  img;

   this(string n, int xs, int ys, int t, int maxv) {
      name   = n;
      xsize  = xs;
      ysize  = ys;
      tuple  = t;
      maxval = maxv;
      img    = new int[xs*ys*t];
   }
   
   Iplane clone() {
      int i;
      auto dup = new Iplane(name,xsize,ysize,tuple,maxval);
      for (i=0; i<(xsize*ysize*tuple);i++) {
         dup.img[i] = this.img[i];
      }
      dup._padx = this._padx;
      dup._pady = this._pady;
      return dup;
   }

   Iplane blankClone() {
      Iplane dup = new Iplane(name,xsize,ysize,tuple,maxval);
      dup._padx = this._padx;
      dup._pady = this._pady;
      return dup;
   }

   Iplane blankClone(int fill) {
      auto dup = this.clone();
      dup.img[] = fill;
      return dup;
   }

   Iplane invert() {
     int n;
     Iplane iout = this.clone();
     for (n=0; n<(xsize*ysize*tuple); n++) {
        iout.img[n] = maxval - iout.img[n];
     }
    return iout;
   }

   Iplane pad(int padVal) {
      int row, col;
      int new_xs = xsize + 2 * padVal;
      int new_ys = ysize + 2 * padVal;
      if (tuple != 1) {
         throw new Exception("Tried to pad a multi-tuple Iplane");
      }
      Iplane padded = new Iplane("grey", new_xs, new_ys, 1, maxval);
      for (row=0; row < ysize; row++) {
         for(col=0; col< xsize; col++) {
            padded.img[(row+padVal-1)*new_xs+(col+padVal-1)] = img[row*xsize+col];
         }
      }
      padded._padx = padVal;
      padded._pady = padVal;
      return padded;
   }

   Iplane unpad(int xp, int yp) {
       _padx = xp;
       _pady = yp;
       return this.unpad();
   }

   Iplane unpad() {
     int row, col;
     int new_xs = xsize - 2 * _padx;
     int new_ys = ysize - 2 * _pady;
     if (tuple != 1) {
        throw new Exception("Tried to unpad a multi-tuple Iplane");
     }
     Iplane unpadded = new Iplane("grey", new_xs, new_ys, 1, maxval);
     for (row=0; row < new_ys; row++) {
        for (col=0; col < new_xs; col++) {
           unpadded.img[row*new_xs+col] = img[(row+_pady)*xsize+(col+_padx)];
        }
     }
     return unpadded;
   }

   Iplane[] explode() {
      int i, n;
      _size = xsize * ysize * tuple;
      Iplane[] output = new Iplane[tuple];
      for (n = 0; n < tuple; n++) {
         output[n] = new Iplane("", xsize, ysize, 1, maxval);
         for (i = 0; i < _size; i++) {
            if (n == (i % tuple)) {
              output[n].img[i/tuple] = img[i];
            }
         }
      }
   return output;
   }

   Iplane applyPixel(int function(int) processPixel) {
      static int flag = 0;
      int i;
      Iplane output = new Iplane(name, xsize, ysize, tuple, maxval);
      for(i=0; i<img.length; i++) {
         output.img[i] = processPixel(img[i]);
         if (flag == 0) {
            flag = 1;
         }
      }
      return output;
   }

   int maxPixel() {
     int max = 0;
     int i;
     for (i=0; i <(xsize*ysize*tuple); i++) {
       if (img[i] > max)
       max = img[i];
     }
     return max;
   }

   Iplane renormalize(int newMaxval) {
      int i;
      int oldMax = 0;
      double scale;
      for (i=0; i<(xsize*ysize); i++) {
         if (img[i] > oldMax) { oldMax = img[i]; }
      }
      scale = cast(double) newMaxval / cast(double) oldMax;
      Iplane output = new Iplane("",xsize,ysize,tuple,newMaxval);
      for (i=0; i<(xsize*ysize); i++) {
         output.img[i] = cast(int) floor(cast(double)img[i]*scale+0.5);
         output.img[i] = (output.img[i] < 0)? 0 : output.img[i];
         output.img[i] = (output.img[i] > newMaxval)? newMaxval : output.img[i];
      }
      return output;
   }

   void setPixel(int x, int y, int fill) {
      if ((x >= xsize) || (y >= ysize) || (x < 0) || (y < 0)) {
         return;
      } 
      img[y*xsize+x] = fill;
   } 
}

Iplane compact(Iplane[] ip) {
   int i;
   auto xs     = ip[0].xsize;
   auto ys     = ip[0].ysize;
   auto mv     = ip[0].maxval;
   auto tuple  = cast(int) ip.length;
   Iplane output = new Iplane("", xs, ys, tuple, mv);
   for (i=0; i<(xs*ys*tuple); i++) {
     output.img[i] = ip[i%tuple].img[i/tuple];
   } 
   return output;
}

struct Fread {
   File f;
   ubyte b;
   bool pushed;

   this(File fl) {
     f = fl;
     pushed = false;
     return;
   }

   ubyte readByte()  {
      ubyte[] b_arry;
      b_arry = new ubyte[1];
      if (pushed) {
        pushed = false;
        return b;
      } else {
        b_arry = f.rawRead(b_arry);
        b = b_arry[0];
        return b;
      }
      throw new Exception("Exited readByte in impossible way.");
    }

    void pushByte(ubyte bt) {
      b = bt;
      pushed = true;
    }

}

void writePAM(string filenam, Iplane i) {
  writeln("*** Cannot yet write PAM format!!!");
}

void writePNM(string filename, Iplane i) {
   // generate exception of tuple != 1
   string magic = "P0";                    // magic number, init to invalid value
   if ((i.tuple == 1) && (i.maxval < 65_536)) {
     magic = "P5";
   } else if ((i.tuple == 3) && (i.maxval < 65_536)) {
     magic = "P6";
   } else if (i.tuple == 1) {
     magic = "P3";     
   } else if (i.tuple == 3) {
     magic = "P4";
   } else {
     writePAM(filename, i);
     return;
   }
   int size = i.xsize * i.ysize * i.tuple;
   int count, n;
   ubyte big, small;
   ubyte [] bbuf;
   auto f = File(filename, "w");
   if (i.maxval < 256) {
      f.writefln("%s\n%d %d\n%d", magic, i.xsize, i.ysize, i.maxval);
      bbuf = new ubyte[size];
      for (n = 0; n < size; n++) {
        bbuf[n] = cast(ubyte)i.img[n];
      }
      f.rawWrite(bbuf);
      f.close();
   } else if (i.maxval < 65536) {
       f.writefln("%s\n%d %d\n%d", magic, i.xsize, i.ysize, i.maxval);
       bbuf = new ubyte[size*2];
       count = 0;
       for (n = 0; n < size; n++) {
           big   = cast(ubyte) (i.img[n] / 256);
           small = cast(ubyte) (i.img[n] % 256);
           bbuf[count] = big;
           count += 1;
           bbuf[count] = small;
           count += 1;
       }
       f.rawWrite(bbuf);
       f.close;
     } else {         // too big for 16-bit, write ascii -- this is untested!!! but so far unneeded
     //throw new Exception("writePNM can't yet do ascii.  Bit depth too great.");
       writeln("Writing ascii ppm file.");
       f.writefln("%s\n%d %d\n%d", magic, i.xsize, i.ysize, i.maxval);  
       count = 0;
       for (n = 0; n < size; n++) {
         f.writef("%d ", i.img[n]);
         count += 1;
         if (count >= i.xsize) {
            f.writeln();
            count = 0;
         }
       }     
   }
return;
}


void seekReturn(Fread f) {
    ReadState state;
    ubyte b;

  state = ReadState.Begin;
  while (state == ReadState.Begin) {
    b = f.readByte();
    if (cast(char) b == '\n') {
      return;
    } else {
    }
  }
  return;
}

int seekIntInHeader(Fread f) {
     ubyte b;
     ReadState state; 
     string output;
     int outi;

   state = ReadState.Begin;
   while (!(state==ReadState.Digit)) {
      b = f.readByte();
      switch (state) {
        case ReadState.Begin, ReadState.Newline : {
           if (isDigit(cast(char) b)) {
              state = ReadState.Digit;
              f.pushByte(b);
              break;
           } else if (cast(char) b == '\n') {  
             state = ReadState.Newline;
             break;  
           } else if (isWhite(cast (char) b)) {
              state = ReadState.Whitespace;
              break;
           } else if ( cast(char) b == '#') {
              state = ReadState.Comment;
              break;
           } else {
              throw new Exception("seekIntInHeader found unexpected character");
           }
        } 
        case ReadState.Whitespace : {
           if (isDigit(cast(char) b)) {
             f.pushByte(b);
             state = ReadState.Digit;
             break;
           } else if ( cast(char) b == '#') {
             state = ReadState.Comment;
             break;
           } else if (isWhite(cast(char) b)) {
             state = ReadState.Whitespace;
             break;
           } else if (cast(char) b == '\n') {
               state = ReadState.Newline;
               break;
           } else {
               state = ReadState.Other;
               break;
           }
        }
        case ReadState.Comment : { 
           if (cast(char) b == '\n') { state = ReadState.Newline; }
           break;
        }
        case ReadState.Digit : {
          output = to!(string)(cast(char) b);
          f.pushByte(b);
          state = ReadState.Digit;
          break;
        } 
        case ReadState.Other : {
           throw new Exception("Improperly formatted header.");
        }
        default: { 
           throw new Exception("Improperly formatted header default.");
        }
   }
  }
   // now whitespace, comments skipped and first digit gotten
   // keep adding digits until something else turns up
   
   state = ReadState.Digit;
   while (state == ReadState.Digit) {
       b = f.readByte();
       if (isDigit(cast(char)b)) {
         output = output ~ (cast(char) b);
       } else {
         state = ReadState.Finished;
         f.pushByte(b);
         if (output.length > 0) {
            outi = to!(int)(strip(output));
            return outi;
         } else {
             throw new Exception("Improperly formatted header. Ghost integer, cannot happen.");
         }
       }
   }
  
  throw new Exception("Impossible error exiting seekIntInHeader");
  return 0;
}

int[3] readPNMheader(Fread f) {
     int[3] answer;
     int xs;
     int ys;
     int max;
     xs  = seekIntInHeader(f);
     ys = seekIntInHeader(f);
     max = seekIntInHeader(f);
     answer[0] = xs;
     answer[1] = ys;
     answer[2] = max;
     return answer;
}

void readPNMbody(File f, Iplane img) {
      ubyte b, s;
      ubyte[] buff;
      int i;
      int maxlen;
    
   if (img.maxval < 256) {
      buff = new ubyte[img.xsize * img.ysize * img.tuple];
      buff = cast(ubyte[]) f.rawRead(buff);
      maxlen = img.xsize * img.ysize * img.tuple;
      maxlen = (maxlen > cast(int) buff.length)? cast(int) buff.length : maxlen;
      for (i=0; i < maxlen; i++) {
         img.img[i] = cast(int) buff[i];
      }
   } else if (img.maxval < 65536) {
      buff = new ubyte[img.xsize * img.ysize * img.tuple * 2];
      buff = cast(ubyte[]) (f.rawRead(buff));
      for (i=0; i < img.xsize*img.ysize*img.tuple; i++) {
         img.img[i] = cast(int) (buff[i*2] * 256 + buff[i*2+1]);
      }
   } else {
      throw new Exception("Bit depth of pnm body > 16 bit.  Giving up.");
   }
}

//*** needs exception handling code
Iplane readPNM(string fileName) {
    int n;
    string pNum;
    int xsize;
    int ysize;
    int tuple;
    int maxval;
    File f;
    Fread fr;
    Iplane img;
  
   // open file for reading, fail if error
   auto magicNum = new ubyte[2];
   auto e        = 0;
   auto w        = 0;
   f             = File(fileName, "r");
   // read magic number, fail if error, choose reader for header
   magicNum = cast(ubyte[])f.rawRead(magicNum);   
   pNum     = cast(string) magicNum;
   // read header, set up reader for pixel values
   fr = Fread(f);
   switch (pNum) {
      case "P1" : { writeln("P1 file.  Not yet supported."); break; }
      case "P2" : { writeln("P2 file.  Not yet supported."); break; }
      case "P3" : { writeln("P3 file.  Not yet supported."); break; }
      case "P4" : { writeln("P4 file.  Not yet supported."); break; }
      case "P5" : {
         auto answer = readPNMheader(fr);
         xsize  = answer[0];
         ysize  = answer[1];
         maxval = answer[2];
         tuple = 1;         // P5 confirms that it is a greyscale image
         img    = new Iplane("img",xsize, ysize, tuple, maxval);
         readPNMbody(f, img);
         writefln("xsize= %d;  ysize= %d; maxval= %d; tuple= %d\n", xsize, ysize, maxval, tuple);
         break;
      }
      case "P6" : {
         auto answer = readPNMheader(fr);
         xsize  = answer[0];
         ysize  = answer[1];
         maxval = answer[2];
         tuple  = 3;
         img    = new Iplane("P6image",xsize, ysize, tuple, maxval);
         readPNMbody(f, img);
         tuple = 3;          // P6 confirms that it is a color image
         writefln("xsize= %d;  ysize= %d; maxval= %d; tuple= %d", xsize, ysize, maxval, tuple);
         break;
      }
      default   : {
         writeln("Invalid file type.  Not a pnm/pam.");
         throw new Exception("Not a pnm/pam.");
         break;
      }
   }
   // make appropriate image, fail if maxval >65535
   // fill pixels, warn if image not finished
   // return image
   f.close();
   return img;
}
