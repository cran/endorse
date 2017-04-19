#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "geodist.h"

static const double DEG2RAD = 0.01745329251994329576923690768;
static const double EARTH_RADIUS_METERS = 6372797.56;

double ArcInRadians(Position from, Position to) {
    double latArc  = (from.lat - to.lat) * DEG2RAD;
    double lonArc = (from.lon - to.lon) * DEG2RAD;
    double latH = sin(latArc * 0.5);
    latH *= latH;
    double lontitudeH = sin(lonArc * 0.5);
    lontitudeH *= lontitudeH;
    double tmp = cos(from.lat*DEG2RAD) * cos(to.lat*DEG2RAD);
    return 2.0 * asin(sqrt(latH + tmp*lontitudeH));
}

double DistanceInMeters(Position from, Position to) {
    return EARTH_RADIUS_METERS*ArcInRadians(from, to);
}
