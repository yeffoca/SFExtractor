import json

data = {
    "infitsList": [],
    "bgfitsList": [],
    "bgSubFitsList": [],
    "seCatsList": [],
    "seSegList": [],
    "header": [
        "NUMBER", "X_IMAGE", "Y_IMAGE", "ID_PARENT",
        "X_WORLD", "Y_WORLD", "THETA_IMAGE", "ELLIPTICITY",
        "A_IMAGE", "B_IMAGE", "PETRO_RADIUS"
    ],
    "seCatsClean": [],
    "targetFile": "files/targets.csv"
}

with open("files/data.json", "w") as f:
    json.dump(data, f, indent=4)
