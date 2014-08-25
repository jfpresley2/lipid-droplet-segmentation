
import std.stdio;
import std.getopt;
import std.c.stdlib;
import medianLib;
import image;

struct Parms {
  string inputImg;
  string outputImg;
  int neighborhood = 120;
  int subsample    = 4;
}

void main(string args[]) {
  Iplane[] imgs;
  Iplane medImg;

  Parms parms;
  getopt(args,
  	 "in", &parms.inputImg,
         "out", &parms.outputImg,
         "neighborhood", &parms.neighborhood,
         "subsample", &parms.subsample);
  if ((parms.inputImg == "")||(parms.outputImg == "")) {
    writeln("Usage: median --in img.pgm --out processedimg.pgm [--neighborhood 120 --subsample 4]  ");
    writeln("Insufficient parameters, program terminating.");
    exit(1);
  }
    
  Iplane inImg  = readPNM(parms.inputImg);
  imgs = inImg.explode();
  foreach(ref Iplane i; imgs) {
     auto filt     = new MedianFilter(i, parms.neighborhood, 50);
     medImg = filt.doBaccor(parms.subsample);
     i = medImg;
  }
  medImg = compact(imgs);
  writePNM(parms.outputImg, medImg);
}