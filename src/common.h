#ifndef _DEF_
#define _DEF_

#ifndef MIN_NUM_POINT_PER_CONTOUR
#define MIN_NUM_POINT_PER_CONTOUR 15
#endif 


#define VIZ //define to enable  visualization utilities 


//********************** dist()
template<typename T>
T dist(T x1, T y1, T z1, T x2, T y2, T z2)
{
    T dx, dy, dz;
    dx = x1 - x2;
    dy = y1 - y2;
    dz = z1 - z2;
    dx *= dx;
    dy *= dy;
    dz *= dz;

    return dx + dy + dz;
}
//******************************************************************************


#endif /*_DEF_*/
