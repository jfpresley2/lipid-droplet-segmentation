
extern (C) int poisson(double mean);
extern (C) double gaussian(double mean, double sigma);
extern (C) void solveMatrix(double *output_x, double*m, 
                            int cols, int rows, double*b);