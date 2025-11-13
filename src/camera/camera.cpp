#include "camera.h"

#include <stdexcept>

Camera::Camera(glm::mat4 view_matrix, float aspect_ratio, float height_angle)
    : view_matrix(view_matrix), aspect_ratio(aspect_ratio), height_angle(height_angle),
      viewplane_height(2.f * std::tan(height_angle / 2.f)), viewplane_width(aspect_ratio * viewplane_height)
{
}

glm::mat4 Camera::getViewMatrix() const
{
    return view_matrix;
}

float Camera::getPlaneWidth() const
{
    return viewplane_width;
}

float Camera::getPlaneHeight() const
{
    return viewplane_height;
}

float Camera::getFocalLength() const
{
    // Optional TODO: implement the getter or make your own design
    throw std::runtime_error("not implemented");
}

float Camera::getAperture() const
{
    // Optional TODO: implement the getter or make your own design
    throw std::runtime_error("not implemented");
}
