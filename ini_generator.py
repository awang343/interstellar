# This script generates .ini configuration files for a rendering engine
# based on JSON scene files located in a specified directory.

import os
import glob
import configparser
import sys

ini_dir = f"inis/"
scene_dir = f"scenefiles/"
textures_dir = f"scenefiles/textures/"
output_dir = f"outputs/"

# Ensure output directories exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(ini_dir, exist_ok=True)

# Default template values
canvas = {"width": "1024", "height": "768"}
feature = {
    "acceleration": "true",
    "parallel": "false",
    "shadows": "true",
    "reflect": "true",
    "texture": "true",
    "bump-mapping": "true",
    "parallax-mapping": "false",
    "super-sample": "false",
    "mipmapping": "true",
}

settings = {
    "texture-filter": "\"bilinear\"",
    "bump-map-filter": "\"nearest\"",
    "bump-scale": "20.0",
    "samples-per-pixel": "1",
    "super-sampler-pattern": "\"grid\"", # options: grid, random, stratified
    "maximum-recursive-depth": "4",
    "only-render-normals": "false",
}

# Loop through all JSON files
for json_file in glob.glob(os.path.join(scene_dir, "*.json")):
    base = os.path.splitext(os.path.basename(json_file))[0]

    # Build paths
    scene_path = json_file
    output_path = os.path.join(output_dir, f"{base}.png")
    ini_path = os.path.join(ini_dir, f"{base}.ini")
    output_mipmaps = os.path.join(output_dir, "mipmaps/")

    # Create config
    config = configparser.ConfigParser()

    config["IO"] = {
        "scene": scene_path,
        "output": output_path,
        "output-mipmaps": output_mipmaps,
        "textures": textures_dir,
    }
    config["Canvas"] = canvas
    config["Feature"] = feature
    config["Settings"] = settings

    # Write ini file
    with open(ini_path, "w") as f:
        config.write(f)

    print(f"Generated: {ini_path}")
