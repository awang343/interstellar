#include "raytracer.h"

#include "intersector.h"
#include "lighting.h"
#include "raytracescene.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <optional>
#include <random>
#include <thread>
#include <utility>
#include <vector>

RayTracer::RayTracer(Config config) : m_config(config)
{
}

// ------------------- Rendering -------------------

void RayTracer::sampler(std::vector<std::pair<float, float>> &samples)
{
    if (!m_config.enableSuperSample)
    {
        samples = {{0.5f, 0.5f}};
        return;
    }

    static const int ss_grid_length = std::sqrt(m_config.samplesPerPixel);

    static thread_local std::mt19937 gen(std::random_device{}());
    static thread_local std::uniform_real_distribution<float> dis(0.f, 1.f);

    switch (m_config.superSamplerPattern)
    {
    case SuperSamplerPattern::Grid:
    {
        static const std::vector<std::pair<float, float>> grid_samples = []
        {
            std::vector<std::pair<float, float>> g(ss_grid_length * ss_grid_length);
            for (int i = 0; i < ss_grid_length; i++)
                for (int j = 0; j < ss_grid_length; j++)
                    g[ss_grid_length * i + j] = {(i + 0.5f) / ss_grid_length, (j + 0.5f) / ss_grid_length};
            return g;
        }();
        samples = grid_samples;
        break;
    }

    case SuperSamplerPattern::Stratified:
        for (int i = 0; i < ss_grid_length; i++)
        {
            for (int j = 0; j < ss_grid_length; j++)
            {
                samples[ss_grid_length * i + j] = {(i + dis(gen)) / ss_grid_length, (j + dis(gen)) / ss_grid_length};
            }
        }
        break;

    case SuperSamplerPattern::Random:
        for (int i = 0; i < m_config.samplesPerPixel; i++)
            samples[i] = {dis(gen), dis(gen)};
        break;
    }
}

void RayTracer::render_pixel(RGBA *imageData, RayTraceScene &scene, const int &i, const int &j,
                             unsigned int &last_ray_id, std::vector<unsigned int> &last_visited,
                             std::vector<std::pair<float, float>> &samples)
{
    static const Camera camera = scene.getCamera();

    static const glm::mat4 view_matrix = camera.getViewMatrix();
    static const glm::mat4 inverse_view = glm::inverse(view_matrix);
    static const glm::vec4 world_camera = inverse_view * glm::vec4(0.f, 0.f, 0.f, 1.f);

    glm::vec3 pixel_color(0.f);
    sampler(samples);

    for (const auto &sample : samples)
    {
        // Position of viewplane pixel in camera space
        const glm::vec4 camera_pixel(
            camera.getPlaneWidth() * ((j + sample.first) / scene_width - 0.5f),
            camera.getPlaneHeight() * ((scene_height - i - sample.second) / scene_height - 0.5f), -1.f, 0.f);

        // Camera to viewplane pixel vector in world space
        const glm::vec3 world_dir = inverse_view * camera_pixel;
        const glm::vec3 norm_dir = glm::normalize(world_dir);

        std::optional<Hit> closest_hit;
        if (m_config.enableAcceleration)
        {
            closest_hit = traverseKDTree(scene, world_camera, norm_dir, true, last_ray_id, last_visited);
        }
        else
        {
            closest_hit = checkIntersection(scene, norm_dir);
        }

        if (!closest_hit)
            continue;

        pixel_color += shadePixel(scene, *closest_hit, world_camera, world_dir, last_ray_id, last_visited);
    }
    pixel_color /= std::size(samples);

    imageData[i * scene_width + j] = {static_cast<uint8_t>(255.f * pixel_color.x),
                                      static_cast<uint8_t>(255.f * pixel_color.y),
                                      static_cast<uint8_t>(255.f * pixel_color.z)};
}

void RayTracer::render(RGBA *imageData, RayTraceScene &scene)
{
    scene_width = scene.getWidth();
    scene_height = scene.getHeight();

    const int max_threads = std::thread::hardware_concurrency() - 2;
    const int num_threads = m_config.enableParallelism ? max_threads : 1;

    if (m_config.enableAcceleration)
    {
        auto tree_build_start = std::chrono::high_resolution_clock::now();
        scene.buildKDTreeSAH();
        auto tree_build_end = std::chrono::high_resolution_clock::now();

        std::cout << "KD Tree build time: " << std::chrono::duration<double>(tree_build_end - tree_build_start).count()
                  << " s\n";
    }

    std::vector<std::vector<unsigned int>> thread_last_visited;
    std::vector<std::vector<std::pair<float, float>>> thread_samples;

    for (int i = 0; i < num_threads; i++)
    {
        thread_last_visited.push_back(std::vector<unsigned int>(scene.getShapes().size(), 0));
        thread_samples.push_back(std::vector<std::pair<float, float>>(m_config.samplesPerPixel));
    }

    const int tile_size = 16;
    const int tiles_x = (scene_width + tile_size - 1) / tile_size;
    const int tiles_y = (scene_height + tile_size - 1) / tile_size;
    const int total_tiles = tiles_x * tiles_y;

    std::atomic<int> next_tile(0);

    const auto render_task = [&](int thread_id)
    {
        unsigned int last_ray_id = 0;
        auto &last_visited = thread_last_visited[thread_id];
        auto &samples = thread_samples[thread_id];

        while (true)
        {
            const int tile_idx = next_tile.fetch_add(1);
            if (tile_idx >= total_tiles)
                break;

            const int ty = (tile_idx / tiles_x) * tile_size;
            const int tx = (tile_idx % tiles_x) * tile_size;

            for (int i = ty; i < std::min(ty + tile_size, scene_height); ++i)
            {
                for (int j = tx; j < std::min(tx + tile_size, scene_width); ++j)
                {
                    render_pixel(imageData, scene, i, j, last_ray_id, last_visited, samples);
                }
            }
        }
    };

    auto render_start = std::chrono::high_resolution_clock::now();

    std::vector<std::thread> workers;
    for (int t = 0; t < num_threads; ++t)
        workers.emplace_back(render_task, t);

    for (auto &w : workers)
        w.join();

    auto render_end = std::chrono::high_resolution_clock::now();
    std::cout << "Render time: " << std::chrono::duration<double>(render_end - render_start).count() << " s\n";
}
