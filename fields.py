import data
import csv
from astropy import units as u
import utilities as utils
import galaxies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
import time

class Field(object):

    # Initializes the field with fits image
    def initWithField(self, infits):
        pass

    # Finds galaxies in field from target list
    def identify(self):
        pass

class RegularField(Field):
    def __init__(self, infits):
        self.infits = infits  # Image of field
        self.targetFile = data.targetFile  # File path for list of targets

    def identify(self):
        self.fieldGalaxies = []  # Dict of galaxies found in this field
        with open(self.targetFile, 'r') as targets:
            reader = csv.DictReader(targets, delimiter=',')
            for row in reader:
                # Pulling RA and dec from targets file
                hdfra = float(row['RA']) * u.deg
                hdfdec = float(row['Dec']) * u.deg

                # Looking for the galaxy in the field
                galaxy, err, bckgrnd = utils.isolateGal(hdfra, hdfdec, self.infits)

                # Ignoring objects not in field and saving objects in field to fieldTargets list
                if galaxy is [0]:
                    continue
                else:
                    self.fieldGalaxies.append(row['No.'])

class SExtractorField(Field):
    def __init__(self, infits, bgfits, seCat, seSeg, seBackSub):
        self.infits = infits  # Image of field
        self.bgfits = bgfits
        self.targetFile = data.targetFile  # File path for list of targets
        self.seCat = seCat
        self.seSeg = seSeg
        self.seBackSub = seBackSub
        self.probChildren = ['167', '181', '190', '142',
                             '144', '184', '204', '207',
                             '220', '150', '152', '156',
                             '161', '194']
        self.galIDList = []

    def identify(self):
        self.fieldGalaxies = {}  # Dict of galaxies found in this field
        with open(self.targetFile, 'r') as targets:
            reader = csv.DictReader(targets, delimiter=',')
            for row in reader:
                # Pulling ID, RA, dec, z, and object name from targets file
                ID = row['No.']
                hdfra = float(row['RA']) * u.deg
                hdfdec = float(row['Dec']) * u.deg
                z = row['Redshift (z)']
                objName = row['Object Name']

                # if ID in self.probChildren:
                #     continue

                # Looking for the galaxy in the field
                galaxy, err, bckgrnd = utils.isolateGal(hdfra, hdfdec, self.infits, self.bgfits)
                galStack = np.dstack((galaxy, err, bckgrnd))
                # plt.imshow(galaxy)
                # plt.show()
                # plt.clf()

                # Ignoring objects not in field and saving objects in field to fieldTargets list
                if galaxy is [0] or len(galaxy) <= 1:
                    continue
                else:
                    self.fieldGalaxies[ID] = (galStack, ID, row['RA'], row['Dec'], z, objName)

    def correlateSE(self):
        self.seIDs = []
        self.seArr = pd.read_csv(self.seCat)
        for galTup in self.fieldGalaxies:
            galStack, ID, RA, dec, z, objName = self.fieldGalaxies[galTup]
            # print(galStack.shape)
            # print(ID)
            galaxy = galaxies.SEGalaxy(ID, RA, dec, z, objName)
            seID, seParentID = galaxy.correlateSE(galStack, self.seArr, self.infits, self.bgfits, iso=True)
            if seID == 0 and seParentID == 0:
                seID, seParentID = galaxy.correlateSE(galStack, self.seArr, self.infits, self.bgfits, iso=False)
            if seID == 0 and seParentID == 0:
                continue

            self.seIDs.append([str(seID), str(seParentID), ID])
            self.galIDList.append(ID)

        self.seIDs = np.array(self.seIDs)

    def mask(self):
        if len(self.seIDs) == 0:
            print('No galaxies found')
        else:
            # Load segmentation image
            seg_map = fits.getdata(self.seSeg)

            # Load SExtractor catalog with safe types
            data = pd.read_csv(self.seCat, header=None, dtype=str)
            seID_df = pd.DataFrame(self.seIDs, columns=[0, 3, 'ID']).astype(str)

            # Merge on OBJECT_ID and PARENT_ID (assumed to be columns 0 and 3)
            filtered_df = pd.merge(data, seID_df, on=[0, 3])
            self.filtered_objects = filtered_df.to_numpy()
            keep_ids = self.filtered_objects[:, 0].astype(int)  # Assuming first column is OBJECT_ID

            # Mask all objects except the selected ones
            mask = np.isin(seg_map, keep_ids)
            seg_map[~mask] = 0  # Set unwanted objects to zero

            # Load original FITS file
            with fits.open(self.seBackSub) as hdul:
                # Make a full deep copy of all HDUs
                hdul_copy = fits.HDUList([hdu.copy() for hdu in hdul])

                original_header = hdul[1].header  # Save header

                # Apply the mask to the data in HDU[1]
                image = hdul_copy[1].data
                self.filtered_image = image * mask  # Set unwanted pixels to 0

                # Create a Primary HDU
                hdu = fits.PrimaryHDU(self.filtered_image)

                # Create an HDU list and write it
                print(f'{self.seCat[-20:-12]}filtered_seg.fits')
                hdulist = fits.HDUList([hdu])
                hdulist.writeto(f'SExtractorStuff/Filtered_Segs/{self.seCat[-20:-12]}filtered_seg.fits', overwrite=True)

    def getAsymmetry(self):
        asymmList = []
        zList = []
        petroMList = []
        for i in range(len(self.filtered_objects)):
            seID, pixRA, pixDEC, seParentID, *_ = self.filtered_objects[i]
            galStack, ID, RA, dec, z, objName = self.fieldGalaxies[self.filtered_objects[i][-1]]
            zList.append(float(z))

            galaxy = utils.isolateGalPix(float(pixRA), float(pixDEC), self.filtered_image, 100, 100)
            galStack = np.dstack((galaxy, galaxy))  # Second parameter should be error
            galaxy = galaxies.SEGalaxy(ID, RA, dec, z, objName)
            absA, err = galaxy.asymmetry(galStack)
            asymmList.append(absA)

            petroMajor, petroMinor = galaxy.getPetro(self.seArr, seID, seParentID)
            petroMList.append(petroMajor)

        return (asymmList, zList, petroMList)

    def showImages(self, id, path, savefig=False):
        for i in range(len(self.filtered_objects)):
            seID, pixRA, pixDEC, *_ = self.filtered_objects[i]
            petroMajor = float(self.filtered_objects[i][-3]) * float(self.filtered_objects[i][-2])
            galStack, ID, RA, dec, z, objName = self.fieldGalaxies[self.filtered_objects[i][-1]]
            galaxy = utils.isolateGalPix(float(pixRA), float(pixDEC), self.filtered_image, 100, 100)

            if id == ID:
                time.sleep(7)
                fig, axes = plt.subplots(1, 2, figsize=(10, 5))  # 1 row, 2 columns

                # First image
                axes[0].imshow(galStack[:, :, 0])
                axes[0].set_title(f'{ID} Original')
                axes[0].axis('off')

                # Second image
                axes[1].imshow(galaxy, vmin=0, vmax=1)
                axes[1].set_title(f'{ID} SE')
                axes[1].axis('off')

                if savefig:
                    plt.savefig(path + '/gal' + ID)
                plt.show()
                plt.clf()

    def listGalaxies(self):
        return self.galIDList




