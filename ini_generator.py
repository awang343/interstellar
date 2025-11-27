import os
import glob
import configparser
import sys

# Directories
name = sys.argv[1]

scene_dir = f"scenefiles/{name}"
output_dir = f"outputs/{name}"
ini_dir = f"inis/{name}"

# Ensure output directories exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(ini_dir, exist_ok=True)

# Default template values
canvas = {"width": "1024", "height": "768"}
feature = {
    "shadows": "true",
    "reflect": "true",
    "refract": "false",
    "texture": "true",
    "texture-filter": "nearest",
    "parallel": "true",
    "super-sample": "false",
    "post-process": "false",
    "acceleration": "true",
    "depthoffield": "false",
}
settings = {"samples-per-pixel": "1", "maximum-recursive-depth": "4"}

# Loop through all JSON files
for json_file in glob.glob(os.path.join(scene_dir, "*.json")):
    base = os.path.splitext(os.path.basename(json_file))[0]

    # Build paths
    scene_path = json_file
    output_path = os.path.join(output_dir, f"{base}.png")
    ini_path = os.path.join(ini_dir, f"{base}.ini")

    # Create config
    config = configparser.ConfigParser()

    config["IO"] = {"scene": scene_path, "output": output_path}
    config["Canvas"] = canvas
    config["Feature"] = feature
    config["Settings"] = settings

    # Write ini file
    with open(ini_path, "w") as f:
        config.write(f)

    print(f"Generated: {ini_path}")
