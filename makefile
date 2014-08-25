
target: project ldspotter ldspotterWater median colocWater

colocWater: colocWater.d image.d watershedLib.d convolve.d circleFit.d wrapper.a wrapper.di spots.d
	dmd image.d circleFit.d convolve.d watershedLib.d spots.d  colocWater.d -L/sw/lib/libgsl.a -L/sw/lib/libgslcblas.a -Lwrapper.a -ofcolocWater -O -inline #-release

median: median.d medianLib.d image.d
	dmd image.d medianLib.d median.d -ofmedian

ldspotter: ldspotter.d image.d convolve.d circleFit.d wrapper.a wrapper.di spots.d 
	dmd image.d circlefit.d convolve.d spots.d  ldspotter.d -L/sw/lib/libgsl.a -L/sw/lib/libgslcblas.a -Lwrapper.a -ofldspotter -O -inline #-release

ldspotterWater: ldspotterWater.d image.d convolve.d circleFit.d wrapper.a wrapper.di spots.d watershedLib.d
	dmd image.d circlefit.d convolve.d spots.d watershedLib.d ldspotterWater.d -L/sw/lib/libgsl.a -L/sw/lib/libgslcblas.a -Lwrapper.a -ofldspotterWater -O -inline #-release

watershed: image.d watershed.d
	dmd image.d watershed.d -ofwatershed

project: image.d medianLib.d project.d
	dmd image.d medianLib.d project.d -ofproject -O -inline -release

gaussRand.a: gaussRand.c
	gcc -Wall  -lgsl -lgslcblas -lm -I/usr/local/include -I/sw/include -c gaussRand.c -o gaussRand.o
	echo "compiled"
	gcc -shared gaussRand.o -lgsl -lgslcblas -lm  -o gaussRand.a

wrapper.a: wrapper.c
	gcc -Wall -lgsl -lgslcblas -lm -I/usr/local/include -I/sw/include -c wrapper.c -o wrapper.o
	echo "wrapper compiled"
	gcc -shared wrapper.o -L/sw/lib -lgsl -lgslcblas -lm -o wrapper.a

#gcc -shared wrapper.o -L/sw/lib/libgslcblas.a -l/sw/lib/libgsl.a -lm -o wrapper.a


