import std.stdio;
import std.c.stdlib;
import std.string;
import std.conv;
import std.getopt;
import std.file;
import image;
import medianLib;
   
string makeName(string root, string extension, int i, int digits) {
   string fmt = "%s%0" ~ std.string.format("%d",digits) ~ "d.%s";
   return std.string.format(fmt, root, i, extension);
}

int findLimit(string root, string extension, int first, int digits) {
   int i = first - 1;
   while (true) {
      string name = makeName(root, extension, ++i, digits);
      //writefln("File found =  %s", name);
      if (!std.file.exists(name)) { 
        i--;
        break; 
      }
   }
   return i;
}

/* 
    Parms has some defaults set that can be overriden from command line.
    These can be edited easily if you prefer different defaults.

 */

struct Parms {          
  string root;
  int first = 0;
  int last = -1;
  int digits = 3;
  string extension;
  string outputfile;  
  int neighborhood = 120; 
  bool nomedian = false;
  int normalize = 4095;
}

void main(string[] args) {
   Parms parms;
   int limit;
   Iplane[] accums;
   getopt(args,
      "root", &parms.root,
      "first", &parms.first,
      "last", &parms.last,
      "digits", &parms.digits,
      "extension", &parms.extension,
      "outputfile", &parms.outputfile,
      "neighborhood", &parms.neighborhood,
      "nomedian", &parms.nomedian,
      "normalize", &parms.normalize,
   );
   if ((parms.root == "") || (parms.outputfile == "")) {
     writeln("Usage: project --root root [--first 0 --last 12 --digits 3 ] --extension ppm --neighborhood 120 --outputfile projected.ppm");
     writeln("No args, program terminating.");
     exit(1);
   }
   string root  = parms.root;
   if (parms.last > 0) {
     limit = parms.last;
   } else {
     limit = findLimit(parms.root, parms.extension, parms.first, parms.digits); 
   }
   if (parms.neighborhood < 0) {
     parms.nomedian = true;
   }
   string ext   = parms.extension;
   string proj  = parms.outputfile;
   bool done    = false;
   int i;
   Iplane accum;
   for (i=0; i<limit; i++) {
     string name = makeName(root, ext, i, parms.digits);
     writeln(name);
     Iplane temp    = readPNM(name);
     Iplane[] tlist = temp.explode();
     if (!parms.nomedian) {
        foreach (ref t; tlist) {
            auto filt   = new MedianFilter(t, parms.neighborhood, 50);
            filt.ssMedian(4);
            t           = filt.doSub();
         }
     }
     temp = compact(tlist);
     if (i==0) {
       accum = temp.clone();
     } else {
       if (   (accum.xsize != temp.xsize)
           || (accum.ysize != temp.ysize) 
           || (accum.tuple != temp.tuple))  {
          throw new Exception("Image sizes don't match.  Projection failed.");
       }
       accum.img[] += temp.img[];
     }

   }
   int maxPixel = accum.maxPixel();
   accums = accum.explode();
   int counter = 0;
   foreach(ref Iplane a; accums) {
     //a         = a.renormalize(4095);
       if (!parms.nomedian) {
         auto filt = new MedianFilter(a, parms.neighborhood, 50);
         filt.ssMedian(4);
	 a         = filt.doSub(); 
       }
       if (parms.normalize < 0) {
         ;
       } else {
         a = a.renormalize(parms.normalize);
       }
   }
   accum = compact(accums);
 
   writePNM(proj, accum);
   writeln(proj);
}