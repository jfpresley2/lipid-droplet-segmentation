import std.stdio;
import std.string;
import std.math;
import image;

Iplane lineConvX(Iplane img, real norm, real[] kern) {
   int row, col, i;
   real sum = 0.0;
   //writePNM("img_before"~img.name~".pgm",img);
   auto xs    = cast(int) img.xsize;
   auto ys    = img.ysize;
   auto  m    = img.maxval;
   auto imgOut = new Iplane(img.name,xs,ys,1,m);
   auto lim1 = kern.length;
   auto lim2 = xs - lim1 - 1;
   auto centr = lim1 / 2;
   for (row = 0; row < (ys - 1); row++) {
     for (col = 0; col < lim2; col++) {
        sum = 0.0;
        for (i = 0; i < kern.length; i++) { 
           sum += kern[i] * cast(real) img.img[row*xs+col+i];
        }
           imgOut.img[row*xs+col+centr] = cast(int) floor(sum / norm + 0.5);
     }
   }   
   //writePNM("img"~img.name~".pgm",imgOut);
   return imgOut;
}

Iplane lineConvY(Iplane img, real norm, real[] kern) {
   int row, col, i;
   real sum   = 0.0;
   auto xs    = img.xsize;
   auto ys    = img.ysize;
   auto  m    = img.maxval;
   auto imgOut = new Iplane(img.name,xs,ys,1,m);
   auto lim1 = kern.length;
   auto lim2 = ys - lim1 - 1;
   auto centr = lim1 / 2;
   for (row = 0; row < lim2; row++) {
     for (col = 0; col < (xs - 1); col++) {
        sum = 0.0;
        for (i = 0; i < kern.length; i++) { 
           sum += kern[i] * cast(real) img.img[(row+i)*xs+col];
        }
        imgOut.img[(row+centr)*xs+col] = cast(int) floor(sum / norm + 0.5);
     }
   }   
   return imgOut;
}

real gauss2(real sigma, real x) {
  return (1.0 / sqrt(2.0 * 3.1415926) * sigma) * exp(-(x*x / (2.0 * sigma*sigma)));
}

real gauss3(real sigma, real x, real y) {
  return (1.0 / 2.0 * 3.1415926) * sigma*sigma * exp(-(x*x + y*y)) / (2.0 * sigma*sigma);
}

struct Kvals {
  real norm;
  real[] kern;
}

auto make2DGauss(real sigma) {
   int i;
   Kvals answer;
   int pixWidth = 2 * (cast(int) ceil(sigma * 3.0)) + 1;
   real[] kern = new real[pixWidth];
   int center = (pixWidth - 1) / 2;
   real centR = cast(real) center;
   real ksum = 100.0;
   for(i=0; i<=center; i++) {
      kern[i] = cast(int) floor(((ksum * gauss2(sigma, (centR - i)) + 0.5)));
   }
   for(i=center+1; i<pixWidth; i++) {
      kern[i] = cast(int) floor(ksum * gauss2(sigma, (i - centR)) + 0.5);   
   }
   answer.kern = kern;
   real norm = 0;
   foreach(e; kern) {
      norm += e;
   }
   //writeln("pixWidth of kernel is ", answer.kern);
   answer.norm = norm;
   //writeln("normalizing factor for kernel is ", answer.norm);
   return answer;
}

Iplane gauss_aux(real sigma, Iplane i) {
   real[] kern;
   real norm;
   //writeln("About to make 2DGauss kernel."); 
   Kvals k = make2DGauss(sigma);
   //*** problematic when sigma is too big relative to image
   //writefln("2D kernel successfully made. k.kern.length=%d",k.kern.length);
   auto i2 = lineConvX(i,k.norm,k.kern);
   //writeln("lineConvX completed.");
   auto i3 = lineConvY(i2,k.norm,k.kern);
   //writeln("lineConvY in function gauss_aux completed.\ngauss_aux completed.");
   return i3;
}

Iplane gauss(Iplane i, double sigma) {
  //writeln("gauss entered.");
  auto result = gauss_aux(sigma,i);      

/*
  Iplane[] j = i.explode();
  //writeln("Image exploded into ", j.length, " images.");
  Iplane[] k = new Iplane[j.length];
  foreach (n, elt; j) {
     //writeln("About to gauss_aux");
     k[n] = gauss_aux(sigma, elt);
     //writeln("gauss_aux completed.");
  }
  auto result = compact(k);
  //writeln("Gaussed images compacted.");
*/
  return result;
}