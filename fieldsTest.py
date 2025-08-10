import fields
import json
import matplotlib.pyplot as plt

with open("files/data.json") as f:
    file_data = json.load(f)

infitsList = file_data["infitsList"]  # List of apertures (fits images)
bgfitsList = file_data["bgfitsList"]  # List of background maps
bgSubFitsList = file_data["bgSubFitsList"]  # List of background subtracted fits
seSegList = file_data["seSegList"]  # List of segmentation maps
seCatsClean = file_data["seCatsClean"]  # List of cleaned catalogs

path = 'files/SavedImages/Tophat3'

IDList = []
absAList = []
zList = []
petroMList = []
for i in range(len(infitsList)):
    seField = fields.SExtractorField(infitsList[i], bgfitsList[i],
                                     seCatsClean[i], seSegList[i], bgSubFitsList[i])
    seField.identify()
    seField.correlateSE()
    seField.mask()
    asymmList, redshiftList, petroList = seField.getAsymmetry()
    ids = seField.listGalaxies()

    for i in range(len(ids)-1):
        if ids[i] not in IDList:
            IDList.append(ids[i])
            absAList.append(asymmList[i])
            zList.append(redshiftList[i])
            petroMList.append(petroList[i])
            if ids[i] in ['167', '181', '190', '142',
                             '144', '184', '204', '207',
                             '220', '150', '152', '156',
                             '161', '194']:
                seField.showImages(ids[i], path, savefig=False)
            # seField.showImages(ids[i])

plt.scatter(zList, absAList)
plt.yscale('log')
plt.minorticks_off()
plt.title('absA vs z')
plt.show()
plt.clf()

plt.scatter(petroMList, absAList)
# plt.yscale('log')
# plt.minorticks_off()
plt.title('absA vs r')
plt.show()
plt.clf()
