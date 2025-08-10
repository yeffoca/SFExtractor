import os
import subprocess
import utilities as utils
import json

# Initializing json with data
subprocess.run(["python", "dataConfigurator.py"], check=True)

# Opening json, so it can be populated
with open("files/data.json") as f:
    file_data = json.load(f)

# Various lists that need to be used or populated
# All lists (other than header) are to be filled with file paths
infits = file_data["infitsList"]  # List of apertures (fits images)
bgfitsList = file_data["bgfitsList"]  # List of background maps
bgSubFitsList = file_data["bgSubFitsList"]  # List of background subtracted fits
seCatsList = file_data["seCatsList"]  # List of SE catalogs
seSegList = file_data["seSegList"]  # List of segmentation maps
header = file_data["header"]  # SE catalogue header
seCatsClean = file_data["seCatsClean"]  # List of cleaned catalogs

# Designating config file
default = False
if default:
    config_file = 'default.sex'  # As given by SExtractor gods. Never changes.
else:
    config_file = 'experiments.sex'  # My dumbass attempt at a config file

# Controls whether SExtractor scan is run
# True produces catalogs
buildCats = True

# Controls whether checkimages are produced
checkImages = False

# Iterating through all apertures
fitsList = os.listdir('files/Apertures')
for fits in fitsList:
    infits.append(f'files/Apertures/{fits}')  # Populating infits list

    identifier = fits[8:18].replace('_', '')
    base_dir = os.path.abspath("SExtractorStuff")
    fits_path = os.path.abspath(f"files/Apertures/{fits}")
    cat_dir = os.path.join(base_dir, "Cats")
    cat_path = os.path.join(cat_dir, f"{identifier}CAT.txt")

    seCatsList.append(cat_path)  # Populating cat list

    # Ensure output directory exists
    os.makedirs(cat_dir, exist_ok=True)

    # Build and run the command
    if buildCats:
        cmd = ["sex", fits_path, "-c", config_file, "-CATALOG_NAME", cat_path]
        try:
            subprocess.run(cmd, check=True, cwd=base_dir)
        except subprocess.CalledProcessError as e:
            print(f"Error running SExtractor on {fits}: {e}")

    # Cleaning the catalogues to be more dataframe/array friendly
    utils.catCleaner(cat_path, header)
    seCatsClean.append(f'{cat_path[:-4]}CLEAN.csv')  # Populating clean cats list

    seg_dir = os.path.join(base_dir, "Segmentations")
    seg_path = os.path.join(seg_dir, f"{identifier}SEG.txt")
    cmdSeg = ["sex", fits_path, "-c", config_file, "-CHECKIMAGE_TYPE", "SEGMENTATION", "-CHECKIMAGE_NAME", seg_path]
    seSegList.append(seg_path)  # Populating segmentation list

    back_dir = os.path.join(base_dir, "Backgrounds")
    back_path = os.path.join(back_dir, f"{identifier}BACK.txt")
    cmdBack = ["sex", fits_path, "-c", config_file, "-CHECKIMAGE_TYPE", "BACKGROUND", "-CHECKIMAGE_NAME", back_path]
    bgfitsList.append(back_path)  # Populating background list

    subBack_dir = os.path.join(base_dir, "-Backgrounds")
    subBack_path = os.path.join(subBack_dir, f"{identifier}-BACK.txt")
    cmdSubBack = ["sex", fits_path, "-c", config_file, "-CHECKIMAGE_TYPE", "BACKGROUND", "-CHECKIMAGE_NAME", subBack_path]
    bgSubFitsList.append(subBack_path)  # Populating subtracted background list

    commands = [
        ("SEGMENTATION", cmdSeg),
        ("BACKGROUND", cmdBack),
        ("OBJECTS", cmdSubBack),
    ]

    errors = []

    # Ensure all output directories exist
    output_paths = [seg_path, back_path, subBack_path]
    for path in output_paths:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    if checkImages:
        for label, cmd in commands:
            try:
                subprocess.run(cmd, check=True, cwd=base_dir)
            except subprocess.CalledProcessError as e:
                errors.append(f"Error running {label} on {fits}: {e}")

        # Optionally print/log all errors at the end
        if errors:
            print("The following errors occurred:")
            for err in errors:
                print(err)


# Saving updated lists back to data.json
with open("files/data.json", "w") as f:
    json.dump(file_data, f, indent=4)
