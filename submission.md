## Project 6: Final Project Gear Up

The project handout can be found [here](https://cs1230.graphics/projects/final/gear-up).

### Test Cases

This project implements the relativistic raytracer that traces through a wormhole from the paper [Visualizing Interstellar’s Wormhole](https://arxiv.org/pdf/1502.03809). The camera is placed so that it always faces towards the center of the wormhole. The wormhole connects two universes, which we call "lower celestial sphere" and "upper celestial sphere." Each celestial sphere is equipped with a panoramic texture image that provides the lighting of the raytracer. 

The wormhole has three parameters: rho, a, M.
- rho: a parameter that controls the width of the wormhole (the radius of the sphere)
- a: a parameter that controls the depth of the wormhole
- M: a parameter that controls the distortion around the wormhole.

We mainly test the functionality of the relativistic raytracer and tweak these parameters. We check if the resulting image matches with the expected outcome from the paper. Unless specified, we adhere to the following setting:
- The pictures have resolution 400 x 300.
- The default values for rho is 1.0.
- The default value for a is 0.05.
- The default value for M is 0.02.

Basic Tests:

1. lower_sphere_blank (the lower sphere is white; the upper sphere is black; renders the position of the wormhole)

<img width="400" height="300" alt="lower_sphere_blank" src="https://github.com/user-attachments/assets/aabb3c09-1eda-4886-b25b-d75c4493e459" />

2. lower_sphere_arrow_lower (the lower sphere is a red arrow pointing right; the upper sphere is black; renders the wormhole)

<img width="400" height="300" alt="lower_sphere_arrow_lower" src="https://github.com/user-attachments/assets/ecdd3ab4-86b1-42d5-b03d-256e3df1ec9f" />


3. lower_sphere_arrow_upper (the lower sphere is white; the upper sphere is a red arrow pointing right but the image warps the direction of the arrow; renders the wormhole)

<img width="400" height="300" alt="lower_sphere_arrow_lower" src="https://github.com/user-attachments/assets/62c0037e-aebe-4c1d-b8de-fffa103d78ed" />

4. lower_sphere_arrow_no_wormhole (renders the lower sphere with the arrow, but without the wormhole; achieved by setting a = rho = 0 and M = 0.01)

<img width="400" height="300" alt="lower_sphere_arrow_no_wormhole" src="https://github.com/user-attachments/assets/bf9ca2ec-5b25-4c37-bd74-f9e8d58081c5" />

5. lower_sphere_distortion_lower (the lower sphere consists of rainbow color stripes; the upper sphere is black)

<img width="400" height="300" alt="lower_sphere_distortion_lower" src="https://github.com/user-attachments/assets/3a4f57e0-1019-427b-a166-60b0b05df037" />

7. lower_sphere_distortion_upper (the lower sphere is black; the upper sphere consists of rainbow color stripes)

<img width="400" height="300" alt="lower_sphere_distortion_upper" src="https://github.com/user-attachments/assets/c4997dc0-95a1-414b-b872-4db19dae9e43" />

Rendered Images:

1. saturn (the lower sphere is an image of galaxies; the upper sphere is an image of saturn)

<img width="400" height="300" alt="saturn" src="https://github.com/user-attachments/assets/2ba66ceb-f86c-430f-82d4-237f0382b203" />

2. saturn_high_res (same picture as 1 but with resolution 1600 x 1200)

<img width="1600" height="1200" alt="saturn_high_res" src="https://github.com/user-attachments/assets/aa442793-8e71-47da-a932-81f0cc58d182" />

3. stars (the upper sphere is an image of galaxies)



4. stars_high_res (same picture as 3 but with resolution 1600 x 1200)

<img width="1600" height="1200" alt="stars_high_res" src="https://github.com/user-attachments/assets/925a4b34-5c47-464c-8c80-33190ef434b5" />


5. saturn_long (saturn, a = 5.0)

<img width="400" height="300" alt="saturn_long" src="https://github.com/user-attachments/assets/f77dcd0f-1239-4309-a4ae-61b0f8f03bb0" />


6. saturn_wide (saturn, M = 0.2)

<img width="400" height="300" alt="saturn_wide" src="https://github.com/user-attachments/assets/6247a24d-86db-4e94-9fe4-b8bf4b3ab989" />


7. saturn_thin (saturn, rho = 0.5)

<img width="400" height="300" alt="saturn_thin" src="https://github.com/user-attachments/assets/5a8f9f34-a4ce-4311-917d-821d73b6c6bd" />


### Design Choices

I stored all relevant info of the rendered scene into json files and wrote a sceneparser helper function to help parse the scene. I also separated the numerical method of solving ODE into a helper class. 

### Collaboration/References

This work implements the paper [Visualizing Interstellar’s Wormhole](https://arxiv.org/pdf/1502.03809). 
I have used chatgpt to write json parser and the QImage generator for me. 

### Known Bugs

For the higher resolution pictures, I need to use smaller time steps for my numerical integration method. However, this would compromise the rendering time. I believe that with enough time I could render very accurate images. 

### Extra Credit
