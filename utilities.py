import os
import math
import datetime
import numpy as np
from astropy.modeling.models import Gaussian2D
# from photutils.datasets import make_noise_image
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import matplotlib.pyplot as plt
import time
import csv


def dms_to_degrees(dms_string):     # FOR USE ON DECLINATION AND LATITUDE
    degrees, minutes, seconds = map(float, dms_string.split(':'))
    degrees += minutes / 60 + seconds / 3600
    return degrees

def time_to_degrees(time_string):   # FOR USE ON RIGHT ASCENSION AND HOUR ANGLE
    hours, minutes, seconds = map(float, time_string.split(':'))
    total_seconds = hours * 3600 + minutes * 60 + seconds
    degrees = total_seconds / 240
    return degrees

def degrees_to_time(degree_string):   # FOR USE ON RIGHT ASCENSION AND HOUR ANGLE
    totalSecs = float(degree_string) * 240
    time = str(datetime.timedelta(seconds=totalSecs))
    return time

def degrees_to_dms(degreeStr):     # FOR USE ON DECLINATION AND LATITUDE
    degrees, decimals = map(int, degreeStr.split('.'))
    decimals = str(decimals * 60)
    minutes = decimals[0:2]
    seconds = str(int(decimals[2:]) * 60)
    return f'{degrees}:{minutes}:{seconds[0:2]}.{seconds[2:]}'

def stringSeparator(str, interval=2, separator=':'):
    return separator.join([str[i:i + interval] for i in range(0, len(str), interval)])

# Shifts an individual wavelength based on a given value for z
def blueshifter(z, wavelength):  # in microns!!!
    restWL = wavelength / (z + 1)
    return restWL

# Shifts a list of wavelengths based on a given value for z
def listBlueshifter(z, wavelengths):  # in microns!!!
    restList = []
    for WL in wavelengths:
        restWL = WL / (z + 1)
        restList.append(restWL)
    return restList

# Avoids hidden folders/files when iterating through a directory
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f

def errorSum(array):
    array *= array
    sum = np.nansum(array)

    return np.sqrt(sum)

# Theta in degrees
# def mockGalaxy(xSD, ySD, theta, nx, ny):
#     g = Gaussian2D(100.0, nx/2, ny/2, xSD, ySD, theta=theta * np.pi / 180.0)
#     y, x = np.mgrid[0:ny, 0:nx]
#     noise = make_noise_image((ny, nx), distribution='gaussian', mean=0.0,
#                              stddev=2.0, seed=1234)
#     data = g(x, y) + noise
#
#     return data

def isolateGal(hdfra, hdfdec, infits, bgfits):
    infile = infits.removesuffix(".fits") + ".fits"
    bgfile = bgfits.removesuffix(".fits") + ".fits"

    # Open the image - the data and the header separately
    fitsimage = fits.open(infile)
    ncards = len(fitsimage)
    headerword = 'CRVAL1'

    phu = 0
    sciii = phu
    badc = 0
    for ii in range(0, ncards):
        headname = fitsimage[ii].name
        try:
            valhead = fitsimage[ii].header[headerword]
            sciii = ii
            break
        except:
            badc += 1
            valhead = "INDEF"

    headersci = fitsimage[sciii].header

    # Now grab the necessary info from the header (usually [SCI] or [1] extension) that is relevant
    # to the world-coordinate-system wcs
    wcs = WCS(headersci)

    xHDF, yHDF = wcs.all_world2pix([[hdfra / u.deg, hdfdec / u.deg]], 0)[0]

    if xHDF < 0 or yHDF < 0 or np.isnan(xHDF) or np.isnan(yHDF):
        return [0], [0], [0]

    xRange = 70
    yRange = 70

    xLower = round(xHDF - xRange)
    xUpper = round(xHDF + xRange-1)
    yLower = round(yHDF - yRange)
    yUpper = round(yHDF + yRange-1)

    with fits.open(infits) as hdul:
        data = hdul[1].data
        err = hdul['ERR'].data
        # err = data.copy()
    with fits.open(bgfits) as bghdul:
        bgdata = bghdul[1].data
        bgerr = bghdul['ERR'].data

    dataYMax, dataXMax = data.shape
    if yUpper < dataYMax and xUpper < dataXMax:
        galaxy = data[yLower:yUpper, xLower:xUpper]
        err = err[yLower:yUpper, xLower:xUpper]
        bckgrnd = bgdata[yLower:yUpper, xLower:xUpper]

        return galaxy, err, bckgrnd
    else:
        return [0], [0], [0]

def isolateGalPix(x, y, image, xRange, yRange):
    xLower = round(x - xRange)
    xUpper = round(x + xRange-1)
    yLower = round(y - yRange)
    yUpper = round(y + yRange-1)

    data = image

    dataYMax, dataXMax = data.shape
    if yUpper < dataYMax and xUpper < dataXMax:
        galaxy = data[yLower:yUpper, xLower:xUpper]
        # err = err[yLower:yUpper, xLower:xUpper]

        return galaxy
    else:
        print("No galaxy found")
        return [0]

def cropGalaxy(galaxy, err, bg, centerIso):
    (yMax, xMax) = galaxy.shape

    if centerIso == True:
        originXY = np.where(galaxy == galaxy[round(yMax/2 - 20):round(yMax/2 + 20),
                                      round(xMax/2 - 20):round(xMax/2 + 20)].max())
    else:
        originXY = np.where(galaxy == galaxy.max())

    center = (originXY[1][0], originXY[0][0])
    dimensions = (yMax, xMax, center)

    # plt.imshow(galaxy[round(yMax/2 - 20):round(yMax/2 + 20),
    #                               round(xMax/2 - 20):round(xMax/2 + 20)])
    # plt.show()
    # plt.clf()



    # Finding distance in pixels between original center and new center
    originalCenter = find_center(galaxy)
    # plt.imshow(galaxy)
    # plt.scatter(originalCenter[0], originalCenter[1], color='red', s=5)
    # plt.title('original')
    # plt.show()
    # plt.clf()
    #
    # plt.imshow(galaxy)
    # plt.scatter(center[0], center[1], color='red', s=5)
    # plt.title('new')
    # plt.show()
    distance = np.sqrt((originalCenter[0] - center[1]) ** 2 + (originalCenter[1] - center[0]) ** 2)
    # print(distance)
    if distance > 20:
        return np.array([]), np.array([]), np.array([]), np.array([])

    origToEdgeX = xMax - center[0]
    origToEdgeY = yMax - center[1]

    # Finding margins around center point
    if origToEdgeX < center[0]:
        xMargin = round(origToEdgeX * 0.9)
    else:
        xMargin = round(center[0] * 0.9)

    if origToEdgeY < center[1]:
        yMargin = round(origToEdgeY * 0.9)
    else:
        yMargin = round(center[1] * 0.9)

    croppedArr = galaxy[center[1] - yMargin:center[1] + (yMargin + 1),
                      center[0] - xMargin:center[0] + (xMargin + 1)]

    croppedErr = err[center[1] - yMargin:center[1] + (yMargin + 1),
                      center[0] - xMargin:center[0] + (xMargin + 1)]

    croppedBG = bg[center[1] - yMargin:center[1] + (yMargin + 1),
                 center[0] - xMargin:center[0] + (xMargin + 1)]

    yMax, xMax = croppedArr.shape
    originXY = np.where(croppedArr == croppedArr[round(yMax/2 - 20):round(yMax/2 + 20),
                                      round(xMax/2 - 20):round(xMax/2 + 20)].max())

    center = (originXY[1][0], originXY[0][0])

    return croppedArr, croppedErr, croppedBG, center

def deg_to_xy(hdfra, hdfdec, infits):
    infile = infits.removesuffix(".fits") + ".fits"

    hdfra = float(hdfra) * u.deg
    hdfdec = float(hdfdec) * u.deg

    # Open the image - the data and the header separately
    fitsimage = fits.open(infile)
    ncards = len(fitsimage)
    headerword = 'CRVAL1'

    phu = 0
    sciii = phu
    badc = 0
    for ii in range(0, ncards):
        headname = fitsimage[ii].name
        try:
            valhead = fitsimage[ii].header[headerword]
            sciii = ii
            break
        except:
            badc += 1
            valhead = "INDEF"

    headersci = fitsimage[sciii].header

    # Now grab the necessary info from the header (usually [SCI] or [1] extension) that is relevant
    # to the world-coordinate-system wcs
    wcs = WCS(headersci)

    xHDF, yHDF = wcs.all_world2pix([[hdfra / u.deg, hdfdec / u.deg]], 0)[0]

    return xHDF, yHDF

def check_image(galImage):
    bool = True
    # plt.imshow(galImage)
    # plt.title('checkimage')
    # plt.show()
    if len(galImage) > 1:
        yMax, xMax = galImage.shape
    else:
        return False

    if yMax == 0 or xMax == 0:
        bool = False
    elif 0.0 in galImage:  # Checking if image includes the border
        bool = False
    elif yMax / xMax > 2 or xMax / yMax > 2:  # Checking if the image is cropped strange
        bool = False

    return bool

def find_center(array_2d):
    rows = len(array_2d)
    cols = len(array_2d[0]) if rows > 0 else 0

    center_row = rows // 2
    center_col = cols // 2

    return (center_row, center_col)

def catCleaner(file, header):
    with open(file, 'r') as reader, open(f'{file[:-4]}CLEAN.csv', 'w') as writer:

        # Removing first 11 lines
        count = 0
        while count < 11:
            reader.readline()
            count += 1

        # writer = csv.writer(writer)
        # writer.writerow(data.header)
        rem_lines = reader.readlines()  # read all the remaining lines
        rem_lines_split = []
        for line in rem_lines:
            rem_lines_split.append(line.split())

        writer = csv.writer(writer)
        writer.writerow(header)
        writer.writerows(rem_lines_split)  # write all the remaining lines to another file
