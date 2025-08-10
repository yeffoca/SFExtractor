# fits files with desired aperture
infitsList = ['files/Apertures/jw02767-o002_t001_nircam_clear-f115w_i2d.fits',
              'files/Apertures/jw03547-o006_t002_nircam_clear-f115w_i2d.fits',
              'files/Apertures/jw05594-o105_t107_nircam_clear-f322w2_i2d.fits',
              'files/Apertures/jw01869-o006_t004_nircam_clear-f277w_i2d.fits']
# infitsList = ['files/Apertures/hst_06337_06_wfpc2_f814w_wf_drz.fits', 'files/Apertures/hst_06337_06_wfpc2_f606w_wf_drz.fits', 'files/Apertures/hst_06337_06_wfpc2_f300w_wf_drz.fits']
# infitsList = []

bgfitsList = ['SExtractorStuff/o002back.fits',
              'SExtractorStuff/o006back.fits',
              'SExtractorStuff/o105back.fits',
              'SExtractorStuff/t004back.fits']

bgSubFitsList = ['SExtractorStuff/o002-BACK.fits',
                 'SExtractorStuff/o006-BACK.fits',
                 'SExtractorStuff/o105-BACK.fits',
                 'SExtractorStuff/t004-BACK.fits']

# txt files with sextractor data for each aperture
seCatsList = ['SExtractorStuff/o002CAT.txt',
              'SExtractorStuff/o006CAT.txt',
              'SExtractorStuff/o105CAT.txt',
              'SExtractorStuff/t004CAT.txt']

seSegList = ['SExtractorStuff/o002SEG.fits',
              'SExtractorStuff/o006SEG.fits',
              'SExtractorStuff/o105SEG.fits',
              'SExtractorStuff/t004SEG.fits']

# SExtractor catologue header
header = ['NUMBER', 'X_IMAGE', 'Y_IMAGE',
          'ID_PARENT', 'X_WORLD', 'Y_WORLD',
          'THETA_IMAGE', 'ELLIPTICITY',
          'A_IMAGE', 'B_IMAGE', 'PETRO_RADIUS']

# Cleaned SExtractor catalogs as csv
seCats = ['SExtractorStuff/o002CATCLEAN.csv',
          'SExtractorStuff/o006CATCLEAN.csv',
          'SExtractorStuff/o105CATCLEAN.csv',
          'SExtractorStuff/t004CATCLEAN.csv']

# File with list of targets and their coords
targetFile = 'files/targets.csv'

# seCatalog = 'SExtractorStuff/catalog.txt'

bads = ['172', '164', '180']

