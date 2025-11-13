#pragma once

#include <glm/glm.hpp>

// A class representing a virtual camera.

// Feel free to make your own design choices for Camera class, the functions
// below are all optional / for your convenience. You can either implement and
// use these getters, or make your own design. If you decide to make your own
// design, feel free to delete these as TAs won't rely on them to grade your
// assignments.

class Camera
{
  public:
    Camera(glm::mat4 view_matrix, float aspect_ratio, float height_angle);

    // Returns the view matrix for the current camera settings.
    // You might also want to define another function that return the inverse of
    // the view matrix.
    glm::mat4 getViewMatrix() const;

    float getPlaneWidth() const;
    float getPlaneHeight() const;

    // Returns the focal length of this camera.
    // This is for the depth of field extra-credit feature only;
    // You can ignore if you are not attempting to implement depth of field.
    float getFocalLength() const;

    // Returns the focal length of this camera.
    // This is for the depth of field extra-credit feature only;
    // You can ignore if you are not attempting to implement depth of field.
    float getAperture() const;

  private:
    glm::mat4 view_matrix;
    glm::vec4 world_camera;
    float aspect_ratio;
    float height_angle;
    float viewplane_height;
    float viewplane_width;
};
