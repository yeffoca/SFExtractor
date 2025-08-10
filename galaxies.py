import utilities as utils
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import scipy
import data
import time
import pandas as pd
import sys

class Galaxy(object):

    # Initializes the galaxy using the ID
    # Pulls coordinates
    def __init__(self, ID, RA, DEC, z, obName):
        pass

    # Finds galaxy from fits file
    # Crops image and error around center
    # Option to display image
    def isolate(self, display=False):
        pass

    # Calculates asymmetry
    # Options for absolute method, squared method,
    # and changing the center to get different A values
    def asymmetry(self, abso=True, sqrd=False, multi=False):
        pass

    # Function for calculating absolute asymmetry method
    def absoluteA(self, croppedArr, arr180, croppedErr, err180):
        pass

    # Function for calculating squared asymmetry method
    def squaredA(self, croppedArr, arr180, croppedErr, err180):
        pass

    # Plots image
    def showImage(self):
        pass

    # Returns RA and DEC
    def getCoords(self):
        pass

    # Returns z(redshift)
    def getRedshift(self):
        pass

    # Returns the ID
    def getID(self):
        pass

    # Returns object name
    def getName(self):
        pass

    # Returns dimensions of the cropped image
    def getShape(self):
        pass



class RegularGalaxy(Galaxy):
    def __init__(self, ID, RA, DEC, z, obName):
        self.ID = ID
        self.infitsList = data.infitsList
        self.RA = RA
        self.DEC = DEC
        self.obName = obName
        self.z = z

    def isolate(self, display=False):
        hdfra = float(self.RA) * u.deg
        hdfdec = float(self.DEC) * u.deg

        count = 0
        for infits in self.infitsList:
            # time.sleep(1)
            galaxy, err, bckgrnd = utils.isolateGal(hdfra, hdfdec, infits)
            if len(galaxy) > 1:
                yMax, xMax = galaxy.shape
                if yMax == 0 or xMax == 0:
                    # print(f"Dimension is 0: (ID-{self.ID}, {yMax}, {xMax})")
                    continue
                elif 0.0 in galaxy:
                    continue
                else:
                    originXY = np.where(galaxy == galaxy.max())
                    center = (originXY[1][0], originXY[0][0])
                    (yMax, xMax) = galaxy.shape
                    dimensions = (yMax, xMax, center)

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

                    self.croppedArr = galaxy[center[1] - yMargin:center[1] + (yMargin + 1),
                                 center[0] - xMargin:center[0] + (xMargin + 1)]

                    self.croppedErr = err[center[1] - yMargin:center[1] + (yMargin + 1),
                                 center[0] - xMargin:center[0] + (xMargin + 1)]

                    originXY = np.where(self.croppedArr == self.croppedArr.max())
                    center = (originXY[1][0], originXY[0][0])
                    (yMax, xMax) = self.croppedArr.shape

                    if yMax/xMax > 2 or xMax/yMax > 2:
                        # print(f"Weird dimensions: (ID-{self.ID}, {yMax}, {xMax})")
                        continue
                    else:
                        self.dimensions = (yMax, xMax, center)
                        self.galWithErr = np.dstack((self.croppedArr, self.croppedErr))

                        if 'c1009' == infits[24:29]:
                            self.bckgrnd = data.bgAList[0]
                        elif 'o001' == infits[24:28]:
                            self.bckgrnd = data.bgAList[1]
                        else:
                            self.bckgrnd = data.bgAList[2]

                        return True
            else:
                count += 1
                # print(f'Out of bounds: ID-{self.ID}')
                continue

        #print(f'Did not appear in {count} files')

    def asymmetry(self, abso=True, sqrd=False, multi=False):
        galaxy = self.galWithErr[:, :, 0]
        err = self.galWithErr[:, :, 1]
        yMax, xMax, center = self.dimensions

        if not multi:
            gal180 = scipy.ndimage.rotate(galaxy, 180)
            err180 = scipy.ndimage.rotate(err, 180)
            if abso and sqrd:
                return self.squaredA(galaxy, gal180, err, err180), self.absoluteA(galaxy, gal180, err, err180)
            elif abso:
                return self.absoluteA(galaxy, gal180, err, err180)
            else:
                return self.squaredA(galaxy, gal180, err, err180)
        else:
            r = 1
            centerList = [center, (center[0] + r, center[1]), (center[0] - r, center[1]),
                          (center[0], center[1] + r), (center[0], center[1] - r),
                          (center[0] + r, center[1] + r), (center[0] - r, center[1] + r),
                          (center[0] - r, center[1] - r), (center[0] + r, center[1] - r)]

            AList = []
            errList = []
            for point in centerList:
                origToEdgeX = xMax - point[0]
                origToEdgeY = yMax - point[1]

                if origToEdgeX < center[0]:
                    xMargin = round(point[0] * 0.9)
                else:
                    xMargin = round(point[0] * 0.9)

                if origToEdgeY < point[1]:
                    yMargin = round(origToEdgeY * 0.9)
                else:
                    yMargin = round(point[1] * 0.9)

                croppedArr = galaxy[point[0] - xMargin:point[0] + (xMargin + 1),
                             point[1] - yMargin:point[1] + (yMargin + 1)]

                croppedErr = err[point[0] - xMargin:point[0] + (xMargin + 1),
                             point[1] - yMargin:point[1] + (yMargin + 1)]

                arr180 = scipy.ndimage.rotate(croppedArr, 180)
                err180 = scipy.ndimage.rotate(croppedArr, 180)

                self.height, self.width = croppedArr.shape

                # plt.imshow(abs(croppedArr-arr180))
                # plt.show()

                if abso:
                    A, errA = self.absoluteA(croppedArr, arr180, croppedErr, err180)
                    # A /= np.sqrt(xMax * yMax)
                    AList.append(A)
                    errList.append(errA)

                elif sqrd:
                    sqrdA, errSqrdA = self.squaredA(croppedArr, arr180, croppedErr, err180)
                    # sqrdA /= np.sqrt(xMax * yMax)
                    AList.append(sqrdA)
                    errList.append(errSqrdA)

            AList.append(self.ID)
            errList.append(self.ID)

            return AList, errList

    def absoluteA(self, croppedArr, arr180, croppedErr, err180):
        # print('err:', croppedErr)
        numerator = np.sum(abs(croppedArr - arr180))
        errNum = utils.errorSum(np.sqrt(croppedErr ** 2 + err180 ** 2))

        denominator = 2 * np.sum(abs(croppedArr))
        errDenom = 2 * utils.errorSum(croppedErr)

        asymmetry = numerator / denominator #- self.bckgrnd
        totalErr = np.sqrt((errNum / numerator) ** 2 + (errDenom / denominator) ** 2) * asymmetry

        return (asymmetry, totalErr)

    def squaredA(self, croppedArr, arr180, croppedErr, err180):
        # print(croppedErr)
        numerator = np.sum((croppedArr - arr180) ** 2)
        sub = croppedArr - arr180
        subErr = np.sqrt(croppedErr ** 2 + err180 ** 2)
        squareErr = sub ** 2 * np.sqrt(2 * (subErr / sub) ** 2)
        errNum = utils.errorSum(squareErr)

        denominator = 2 * np.sum(croppedArr ** 2)
        errDenom = 2 * utils.errorSum(np.sqrt(2 * (croppedErr / croppedArr) ** 2) * croppedArr ** 2)

        asymmetry = numerator / denominator #- self.bckgrnd
        # print(asymmetry)
        totalErr = np.sqrt((errNum / numerator) ** 2 + (errDenom / denominator) ** 2) * asymmetry

        return (asymmetry, totalErr)

    def radialAsymm(self, radiiArr, abs=True, sqrd=False):
        self.radiiArr = radiiArr
        gal = self.croppedArr
        galErr = self.croppedErr
        # gal180 = scipy.ndimage.rotate(gal, 180)
        yMax, xMax, center = self.dimensions

        self.asymmList = []
        self.radList = []
        self.errList = []
        for radius in self.radiiArr:
            if radius + center[0] < yMax and radius + center[1] < xMax:
                croppedGal = gal[center[1] - radius:center[1] + radius,
                             center[0] - radius:center[0] + radius]
                croppedErr = galErr[center[1] - radius:center[1] + radius,
                             center[0] - radius:center[0] + radius]
                croppedGal180 = scipy.ndimage.rotate(croppedGal, 180)
                croppedErr180 = scipy.ndimage.rotate(croppedErr, 180)
                if abs:
                    asymm, err = self.absoluteA(croppedGal, croppedGal180, croppedErr, croppedErr180)
                elif sqrd:
                    asymm, err = self.squaredA(croppedGal, croppedGal180, croppedErr, croppedErr180)

                self.asymmList.append(asymm)
                self.radList.append(radius)
                self.errList.append(err)

        # return asymmList, errList

    def plotRadAsymm(self):
        plt.scatter(self.radList, self.asymmList)
        plt.ylabel('Asymmetry')
        plt.xlabel('Radius [pixels]')
        plt.show()
        plt.clf()

        self.showImages()

    def showImages(self, galaxy=True, rotated=False, sub=False):
        arr180 = scipy.ndimage.rotate(self.croppedArr, 180)

        if galaxy:
            time.sleep(5)

            yMax, xMax = self.croppedArr.shape
            plt.imshow(self.croppedArr, origin='lower')#, vmin=0, vmax=1)
            plt.scatter(xMax/2, yMax/2, color='red')
            plt.title(f'{self.ID} Original')
            # plt.savefig(f'{self.ID}_Original.png')
            plt.show()
            plt.clf()

            plt.imshow(self.croppedArr, origin='lower', vmin=0, vmax=0.3)
            plt.scatter(xMax / 2, yMax / 2, color='red')
            plt.title(f'{self.ID} Bright')
            # plt.savefig(f'{self.ID}_Bright.png')
            plt.show()
            plt.clf()

        if rotated:
            plt.imshow(arr180, origin='lower')
            plt.title(f'{self.ID} Rotated')
            plt.show()
            plt.clf()

        if sub:
            plt.imshow(abs(self.croppedArr - arr180) / (2*abs(self.croppedArr)), origin='lower')
            plt.title(f'{self.ID} Subtracted')
            plt.show()
            plt.clf()

    def getCoords(self):
        return self.RA, self.DEC

    def getRedshift(self):
        return self.z

    def getID(self):
        return self.ID

    def getName(self):
        return self.obName

    def getShape(self):
        return self.height, self.width


# Use when utilizing SExtractor data (i.e. center, ellipticity)
class SEGalaxy(Galaxy):
    def __init__(self, ID, RA, DEC, z, obName):
        self.ID = ID
        # self.infitsList = data.infitsList
        # self.bgfitsList = data.bgfitsList
        # self.seCatsList = data.seCats
        # self.seSegList = data.seSegList
        self.RA = float(RA)
        self.DEC = float(DEC)
        self.obName = obName
        self.z = z

    def correlateSE(self, galStack, seArr, infits, bgfits, iso):
        self.seArr = seArr
        galaxy = galStack[:, :, 0]
        err = galStack[:, :, 1]
        bg = galStack[:, :, 2]
        # print(self.ID)
        if np.isnan(galaxy).any():
            nan_indices = np.where(np.isnan(galaxy))
            nanY = nan_indices[0][0]
            nanX = nan_indices[1][0]
            galaxy[nanY, nanX] = bg[nanY, nanX]

        # if self.ID == '158':
        #     plt.imshow(galaxy)
        #     plt.title(f'{self.ID}')
        #     plt.show()

        yMax, xMax = galaxy.shape
        originXY = np.where(galaxy == galaxy[round(yMax/2 - 20):round(yMax/2 + 20),
                                  round(xMax/2 - 20):round(xMax/2 + 20)].max())
        center = (originXY[1][0], originXY[0][0])

        epsilon = 1e-3
        # print(seArr.X_WORLD.head())
        # print(seArr.Y_WORLD.head())
        # print(seArr)
        closeTargets = seArr.loc[(abs(seArr.X_WORLD - self.RA) < epsilon) &
                                 (abs(seArr.Y_WORLD - self.DEC) < epsilon)]

        if not closeTargets.empty:  # Confirm there was at least one close target
            closestTargetIDs = (closeTargets.iloc[0]['NUMBER'], closeTargets.iloc[0]['ID_PARENT'])  # Pulling identifiers of lowest values
            closestTarget = seArr.loc[(seArr.NUMBER == closestTargetIDs[0]) &  # Pulling all information for closest sextractor target
                                      (seArr.ID_PARENT == closestTargetIDs[1])]

            for i in range(len(closeTargets) - 1):
                seRA = list(closeTargets.X_WORLD)[i] * u.deg
                seDec = list(closeTargets.Y_WORLD)[i] * u.deg
                seID = list(closeTargets.NUMBER)[i]
                seParentID = list(closeTargets.ID_PARENT)[i]
                seGal, seErr, seBG = utils.isolateGal(seRA, seDec, infits, bgfits)

                # if self.ID == '167':
                #     print(seID, seParentID)
                #     plt.imshow(seGal, vmin=0, vmax=1)
                #     plt.title('input')
                #     plt.show()
                #     plt.clf()
                # plt.imshow(seGal, vmin=0, vmax=1)
                # plt.title(seID)
                # plt.show()

                if np.isnan(seGal).any():
                    nan_indices = np.where(np.isnan(seGal))
                    nanY = nan_indices[0][0]
                    nanX = nan_indices[1][0]
                    seGal[nanY, nanX] = seBG[nanY, nanX]

                if np.isnan(seGal).any():
                    continue

                seGalCropped, seErrCropped, bgCropped, seCenter = utils.cropGalaxy(seGal, seErr, seBG, iso)
                if seGalCropped is np.array([]):
                    continue

                # print(self.ID)

                #
                # plt.imshow(galaxy)
                # plt.title(f'Galaxy {self.ID}')
                # plt.show()

                #print(self.ID, utils.check_image(seGalCropped))
                closestDiff = 100
                if utils.check_image(seGalCropped):
                    yDiff = abs(seGalCropped.shape[0] - galaxy.shape[0])
                    xDiff = abs(seGalCropped.shape[1] - galaxy.shape[1])
                    sizeDiff = (yDiff + xDiff) / 2

                    if sizeDiff < closestDiff:
                        closestDiff = sizeDiff

                        # print(seGalCropped.shape, seCenter)
                        # print(galaxy.shape, center)
                        # print(seGalCropped.shape)
                        # print(seCenter, center)

                        seYMax, seXMax = seGalCropped.shape
                        smallestYMargin = min(seYMax - seCenter[1], yMax - center[1]) - 1
                        smallestXMargin = min(seXMax - seCenter[0], xMax - center[0]) - 1
                        # print(smallestYMargin, smallestXMargin)

                        if smallestYMargin > seCenter[1] or smallestYMargin > center[1]:
                            continue
                        if smallestXMargin > seCenter[0] or smallestXMargin > center[0]:
                            continue

                        seGalCropped = seGalCropped[seCenter[1] - smallestYMargin:seCenter[1] + smallestYMargin,
                                       seCenter[0] - smallestXMargin:seCenter[0] + smallestXMargin]
                        seErrCropped = seErrCropped[seCenter[1] - smallestYMargin:seCenter[1] + smallestYMargin,
                                       seCenter[0] - smallestXMargin:seCenter[0] + smallestXMargin]
                        # print(center)
                        # print(galaxy.shape)
                        galCropped = galaxy[center[1] - smallestYMargin:center[1] + smallestYMargin,
                                       center[0] - smallestXMargin:center[0] + smallestXMargin]
                        # print(galCropped.shape)
                        errCropped = err[center[1] - smallestYMargin:center[1] + smallestYMargin,
                                     center[0] - smallestXMargin:center[0] + smallestXMargin]

                        # print(seGalCropped.shape)
                        # print(galCropped.shape)
                        # if self.ID == '136':
                        # time.sleep(7)
                        # plt.imshow(seGalCropped)
                        # plt.title(f'SE {self.ID}')
                        # plt.show()
                        # plt.clf()
                        #
                        # plt.imshow(galCropped)
                        # plt.title(f'Original {self.ID}')
                        # plt.show()
                        # plt.clf()
                        asymmetrySE, seTotalErr = self.absoluteA(galCropped, seGalCropped, errCropped, seErrCropped, correlation=True)
                        if asymmetrySE == 0:
                            # time.sleep(7)
                            # plt.imshow(seGalCropped)
                            # plt.title(f'SE {self.ID}')
                            # plt.show()
                            # plt.clf()
                            #
                            # plt.imshow(galCropped)
                            # plt.title(f'Original {self.ID}')
                            # plt.show()
                            # plt.clf()
                            seGalCropped = seGalCropped
                            seID2 = seID
                            seParentID2 = seParentID
                            bckgrnd = bgCropped

                            # plt.imshow(seGalCropped, vmin=0, vmax=1)
                            # plt.title(f'{seID} match')
                            # plt.show()
                            # print(self.ID)
                            # print(seID2, seParentID2)


                            return (seID2, seParentID2)
                            # return True
        return (0, 0)

    def asymmetry(self, galWithErr, abso=True, sqrd=False, multi=False):
        galaxy = galWithErr[:, :, 0]
        err = galWithErr[:, :, 1]
        yMax, xMax = galaxy.shape
        originXY = np.where(galaxy == galaxy.max())
        center = (originXY[1][0], originXY[0][0])


        if not multi:
            gal180 = scipy.ndimage.rotate(galaxy, 180)
            err180 = scipy.ndimage.rotate(err, 180)
            if abso and sqrd:
                return self.squaredA(galaxy, gal180, err, err180), self.absoluteA(galaxy, gal180, err, err180)
            elif abso:
                return self.absoluteA(galaxy, gal180, err, err180)
            else:
                return self.squaredA(galaxy, gal180, err, err180)
        else:
            r = 1
            centerList = [center, (center[0] + r, center[1]), (center[0] - r, center[1]),
                          (center[0], center[1] + r), (center[0], center[1] - r),
                          (center[0] + r, center[1] + r), (center[0] - r, center[1] + r),
                          (center[0] - r, center[1] - r), (center[0] + r, center[1] - r)]

            AList = []
            errList = []
            for point in centerList:
                origToEdgeX = xMax - point[0]
                origToEdgeY = yMax - point[1]

                if origToEdgeX < center[0]:
                    xMargin = round(point[0] * 0.9)
                else:
                    xMargin = round(point[0] * 0.9)

                if origToEdgeY < point[1]:
                    yMargin = round(origToEdgeY * 0.9)
                else:
                    yMargin = round(point[1] * 0.9)

                croppedArr = galaxy[point[0] - xMargin:point[0] + (xMargin + 1),
                             point[1] - yMargin:point[1] + (yMargin + 1)]

                croppedErr = err[point[0] - xMargin:point[0] + (xMargin + 1),
                             point[1] - yMargin:point[1] + (yMargin + 1)]

                arr180 = scipy.ndimage.rotate(croppedArr, 180)
                err180 = scipy.ndimage.rotate(croppedArr, 180)

                self.height, self.width = croppedArr.shape

                # plt.imshow(abs(croppedArr-arr180))
                # plt.show()
                # print(self.bckgrnd)
                if abso:
                    A, errA = self.absoluteA(croppedArr, arr180, croppedErr, err180)
                    # A /= np.sqrt(xMax * yMax)
                    AList.append(A)
                    errList.append(errA)

                elif sqrd:
                    sqrdA, errSqrdA = self.squaredA(croppedArr, arr180, croppedErr, err180)
                    # sqrdA /= np.sqrt(xMax * yMax)
                    AList.append(sqrdA)
                    errList.append(errSqrdA)

            AList.append(self.ID)
            errList.append(self.ID)

            return AList, errList

    def absoluteA(self, croppedArr, arr180, croppedErr, err180, correlation=False):
        self.bg = True
        numerator = np.sum(abs(croppedArr - arr180))
        errNum = utils.errorSum(np.sqrt(croppedErr ** 2 + err180 ** 2))

        denominator = 2 * np.sum(abs(croppedArr))
        errDenom = 2 * utils.errorSum(croppedErr)
        asymmetry = numerator / denominator
        totalErr = np.sqrt((errNum / numerator) ** 2 + (errDenom / denominator) ** 2) * asymmetry

        # if self.bg and not correlation:
        #     bg180 = scipy.ndimage.rotate(bckgrnd, 180)
        #     bgNum = np.sum(abs(bckgrnd - bg180))
        #     asymmetry = (numerator / denominator) - (bgNum / denominator)

        return (asymmetry, totalErr)

    def squaredA(self, croppedArr, arr180, croppedErr, err180, correlation=False):
        self.bg = True
        # print(croppedErr)
        numerator = np.sum((croppedArr - arr180) ** 2)
        sub = croppedArr - arr180
        subErr = np.sqrt(croppedErr ** 2 + err180 ** 2)
        squareErr = sub ** 2 * np.sqrt(2 * (subErr / sub) ** 2)
        errNum = utils.errorSum(squareErr)

        denominator = 2 * np.sum(croppedArr ** 2)
        errDenom = 2 * utils.errorSum(np.sqrt(2 * (croppedErr / croppedArr) ** 2) * croppedArr ** 2)

        asymmetry = numerator / denominator # - self.bckgrnd
        # print(asymmetry)
        totalErr = np.sqrt((errNum / numerator) ** 2 + (errDenom / denominator) ** 2) * asymmetry

        # if self.bg and not correlation:
        #     bg180 = scipy.ndimage.rotate(bckgrnd, 180)
        #     bgNum = np.sum((bckgrnd - bg180) ** 2)
        #     # denominator = np.sum(croppedArr ** 2)
        #     asymmetry = (numerator / denominator) - (bgNum / denominator)
            # print(asymmetry)

        return (asymmetry, totalErr)

    def radialAsymm(self, radiiArr, abs=True, sqrd=False):
        self.radiiArr = radiiArr
        gal = self.croppedArr
        galErr = self.croppedErr
        # gal180 = scipy.ndimage.rotate(gal, 180)
        yMax, xMax, center = self.dimensions

        self.asymmList = []
        self.radList = []
        self.errList = []
        for radius in self.radiiArr:
            if radius + center[0] < yMax and radius + center[1] < xMax:
                croppedGal = gal[center[1] - radius:center[1] + radius,
                             center[0] - radius:center[0] + radius]
                croppedErr = galErr[center[1] - radius:center[1] + radius,
                             center[0] - radius:center[0] + radius]
                croppedGal180 = scipy.ndimage.rotate(croppedGal, 180)
                croppedErr180 = scipy.ndimage.rotate(croppedErr, 180)
                if abs:
                    asymm, err = self.absoluteA(croppedGal, croppedGal180, croppedErr, croppedErr180)
                elif sqrd:
                    asymm, err = self.squaredA(croppedGal, croppedGal180, croppedErr, croppedErr180)

                self.asymmList.append(asymm)
                self.radList.append(radius)
                self.errList.append(err)

        # return asymmList, errList

    def plotRadAsymm(self):
        plt.scatter(self.radList, self.asymmList)
        plt.ylabel('Asymmetry')
        plt.xlabel('Radius [pixels]')
        plt.show()
        plt.clf()

        self.showImages()

    def showImages(self, galaxy=True, rotated=False, sub=False):
        self.croppedGal = self.galWithErr[:, :, 0]
        arr180 = scipy.ndimage.rotate(self.croppedGal, 180)

        if galaxy:
            time.sleep(5)

            yMax, xMax = self.croppedGal.shape
            plt.imshow(self.croppedGal, origin='lower')#, vmin=0, vmax=1)
            plt.scatter(xMax/2, yMax/2, color='red', s=5)
            plt.title(f'{self.ID} Original')
            # plt.savefig(f'{self.ID}_Original.png')
            plt.show()
            plt.clf()

            plt.imshow(self.croppedGal, origin='lower', vmin=0, vmax=1)
            plt.scatter(xMax / 2, yMax / 2, color='red', s=5)
            plt.title(f'{self.ID} Bright')
            # plt.savefig(f'{self.ID}_Bright.png')
            plt.show()
            plt.clf()

            # plt.imshow(self.croppedGal - self.bg, origin='lower', vmin=0, vmax=1)
            # plt.scatter(xMax / 2, yMax / 2, color='red', s=5)
            # plt.title(f'{self.ID} BG Subtracted')
            # # plt.savefig(f'{self.ID}_Bright.png')
            # plt.show()
            # plt.clf()

            if self.seBool:
                yMax, xMax = self.seGalCropped.shape
                plt.imshow(self.seGalCropped, origin='lower', vmin=0, vmax=1)
                plt.scatter(xMax / 2, yMax / 2, color='red', s=5)
                plt.title(f'{self.seID}, {self.seCat} SE {self.seAsymm}')
                # plt.savefig(f'{self.ID}_Bright.png')
                plt.show()
                plt.clf()
            #
            #     yMax, xMax = self.bckgrnd.shape
            #     plt.imshow(self.bckgrnd, origin='lower')#, vmin=0, vmax=0.3)
            #     plt.title(f'{self.ID} SE {self.seAsymm}')
            #     # plt.savefig(f'{self.ID}_Bright.png')
            #     plt.show()
            #     plt.clf()



        if rotated:
            plt.imshow(arr180, origin='lower')
            plt.title(f'{self.ID} Rotated')
            plt.show()
            plt.clf()

        if sub:
            arr180 = abs(self.croppedArr - arr180) #/ (2 * abs(self.croppedArr)))
            # plt.imshow(abs(self.croppedArr - arr180) / (2*abs(self.croppedArr)), origin='lower', vmin=0, vmax=1)
            plt.imshow(arr180, origin='lower', vmin=0, vmax=1)
            plt.title(f'{self.ID} Subtracted')
            plt.show()
            plt.clf()

            print(f'{self.ID}: {np.sum(arr180)}, {np.sum(abs(self.croppedArr))}')

    def saveImagesToFits(self):
        croppedGal = self.galWithErr[:, :, 0]
        gal180 = scipy.ndimage.rotate(croppedGal, 180)

        hduSkySub = fits.PrimaryHDU(croppedGal-self.bg)
        hdulist = fits.HDUList([hduSkySub])
        hdulist.writeto('skysub.fits', overwrite=True)
        hdulist.close()

        hduSkySubAbs = fits.PrimaryHDU(abs(croppedGal - self.bg))
        hdulist = fits.HDUList([hduSkySubAbs])
        hdulist.writeto('skysubabs.fits', overwrite=True)
        hdulist.close()

        hduSub180 = fits.PrimaryHDU(croppedGal - gal180)
        hdulist = fits.HDUList([hduSub180])
        hdulist.writeto('sub180.fits', overwrite=True)
        hdulist.close()

        hduSub180Abs = fits.PrimaryHDU(abs(croppedGal - gal180))
        hdulist = fits.HDUList([hduSub180Abs])
        hdulist.writeto('sub180abs.fits', overwrite=True)
        hdulist.close()

    def getCoords(self):
        return self.RA, self.DEC

    def getRedshift(self):
        return self.z

    def getID(self):
        return self.ID

    def getName(self):
        return self.obName

    def getShape(self):
        return self.height, self.width

    def getPetro(self, seArr, seID, seParentID):
        targetRow = seArr.loc[(seArr.NUMBER == int(seID)) &
                              (seArr.ID_PARENT == int(seParentID))]
        petroRMajor = float(targetRow.PETRO_RADIUS) * float(targetRow.A_IMAGE)
        petroRMinor = float(targetRow.PETRO_RADIUS) * float(targetRow.B_IMAGE)

        return petroRMajor, petroRMinor
