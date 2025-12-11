# Interstellar
Inspired by the wormhole from the movie _Interstellar_ and building off of the ray tracers that we built in CSCI1230, Team Interstellar presents our relativistic ray tracer. 

## How to use it
Within the folder _scenefiles_ there are a number of .jsons that describe how the scene is set up, how the scene moves, and how the wormhole is set up. Simply create a new .json or use an existing one as the command line argument for the program and the program will spit out an image (if only one image is requested) or a folder containing all images (if multiple images are requested).

### Json files fields
| Field | Type| Required | Description|
| --- | --- | --- | --- |
| `upperTexture` | string |Yes | Upper celestial sphere texture |
| `lowerTexture` | string |Yes | Lower celestial sphere texture |
| `outputImage` | string |Yes* | Used when numPhotos = 1|
| `outputFolder`| string| No | Used when numPhotos > 1|
| `rho` |number| Yes*  |Shape parameter of the wormhole|
| `a` |number| Yes* | Spin/throat twist parameter |
| `M` |number| Yes* | Mass parameter (akin to gravity)|
| `outWidth`| int | Yes | Output image resolution|
| `inWidth`| int | Yes | OUtput image resolution|
| viewPlaneWidthAngle| number| Yes | Fov (angle) |
| dt | number | Yes | Integration step|
| cameraDistance | number | Yes* | Default camera position |
| numPhotos| int | Yes | Number of frames|
| cameraPoints | array | no | Camera keyframes, where each entry is `[r, theta, phi, t]`|
| wormholePoints | array | no | Wormhole keyframes, where each entry is `[rho, a, m, t]` |
| objects | array | no | Animated scene objects, where each is an object containing `type`, `size`, `keyframes`|

_* indicates that these are technically not required and may be ignored under certain conditions, but due to a rushed finish and deprecated code, the program will ask that you fill them in_

## Known Bugs

There appear to be some intersection bugs when the cube appears behind the wormhole. `cube_around.json` demonstrates this phenomenon. Additionally, there seems to be some unexpected reflection when a sphere passes by a wormhole. 
