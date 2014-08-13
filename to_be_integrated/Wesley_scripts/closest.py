from __future__ import division
import numpy as np
import numexpr as ne


lonlat = np.array([adcpData.lon[0], adcpData.lat[0]]).T
point_list = np.array([test.lon, test.lat]).T

def closest_point1(points, lon, lat):

    point_list = np.array([lon,lat]).T

    #print point_list
    #print points

    point_list0 = point_list[:, 0]
    points0 = points[:, 0, None]
    point_list1 = point_list[:, 1]
    points1 = points[:, 1, None]
    closest_dist = ne.evaluate('((point_list0 - points0)**2 + \
                    (point_list1 - points1)**2)')

    #closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
    #                (point_list[:, 1] - points[:, 1, None])**2)

    print closest_dist
    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes


def closest_point2(points, lon, lat):

    point_list = np.array([lon,lat]).T


    closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                    (point_list[:, 1] - points[:, 1, None])**2)

    np.diff(abs(point_list[:, 0]), abs(points[:, 0]))

    print closest_dist

    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes

def closest_point_old(points, lon, lat):

    point_list = np.array([lon,lat]).T


    closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                    (point_list[:, 1] - points[:, 1, None])**2)

    print closest_dist

    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes

# WES_COMMENT: FASTEST at 902 microseconds
def closest_point(points, point_list):

    point_list0 = point_list[:, 0]
    points0 = points[:, 0, None]
    point_list1 = point_list[:, 1]
    points1 = points[:, 1, None]

    closest_dist = ((point_list0 - points0) *
                    (point_list0 - points0) +
                    (point_list1 - points1) *
                    (point_list1 - points1)
                    )

    print closest_dist

    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes

