#include <iostream>
#include <cairo/cairo.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
using namespace std;

#define widthOfScreen 400
#define heightOfScreen 400
////////////x-y=0////////////////////
#define lineParameterA 1
#define lineParameterB -1
#define lineParameterC 0
////////////////////////////////
double getRandomDouble(double fmin,double fmax);
int main()
{
    cairo_surface_t *surface;
    cairo_t *cr;
    surface = cairo_image_surface_create (CAIRO_FORMAT_RGB24,widthOfScreen,heightOfScreen);
    cr = cairo_create (surface);
    cairo_new_path (cr);
    cairo_rectangle(cr, 0, 0, widthOfScreen, heightOfScreen);
    cairo_set_source_rgb (cr, 1, 1, 1);
    //cairo_fill(cr);
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,400,400);
    cairo_stroke(cr);
    std::srand(std::time(0));
    std::ofstream ofile("trainingData.txt");
    for(int i=0; i<3000; i++)
    {

        double dCoorX=getRandomDouble(0,widthOfScreen);
        double dCoorY=getRandomDouble(0,heightOfScreen);
        ofile<<"Input: "<<dCoorX<<"  "<<dCoorY;
        ofile<<std::endl;
        cairo_move_to(cr,dCoorX,dCoorY);
        double dTemp=lineParameterA*dCoorX+lineParameterB*dCoorY+lineParameterC;
        bool bFlag=false;
        if(dTemp>=0)
        {
            cairo_set_source_rgb(cr,255,0,0);

            cairo_arc (cr, dCoorX, dCoorY, 2, 0, 2 * M_PI);
            //cairo_fill(cr);
        }
        else
        {
            cairo_set_source_rgb(cr,0,255,0);
            cairo_rectangle(cr,dCoorX,dCoorY,2,2);
            bFlag=true;
            //cairo_fill(cr);

        }
        double iOutPut;
        iOutPut=bFlag==false?0:1;
        ofile<<"Output: "<<iOutPut<<std::endl;
        cairo_stroke(cr);

    }

    ofile.close();
    cairo_set_source_rgb (cr, 1, 1, 1);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

    cairo_surface_write_to_png(surface, "image.png");
    cout << "Hello world!" << endl;
    return 0;
}

double getRandomDouble(double fmin,double fmax)
{
    int iRand=std::rand();
    double temp=fmin+(double)iRand/(double)RAND_MAX*(fmax-fmin);
    return temp;
}
