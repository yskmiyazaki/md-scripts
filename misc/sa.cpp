#include <opencv2/opencv.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

using namespace cv;
using namespace std;

// Function to compute Mean Squared Error between two images
double computeMSE(const Mat& img1, const Mat& img2) {
    Mat diff;
    absdiff(img1, img2, diff);
    diff = diff.mul(diff);
    Scalar sum_diff = sum(diff);
    return (sum_diff[0] + sum_diff[1] + sum_diff[2]) / (img1.total());
}

// Simulated Annealing for optimizing the positions
vector<Point> simulatedAnnealing(const vector<Mat>& images, const vector<Point>& initial_positions, const int width, const int height) {
    vector<Point> current_positions = initial_positions;
    vector<Point> best_positions = initial_positions;
    double best_overlap_error = numeric_limits<double>::max();
    
    // Parameters for simulated annealing
    double temperature = 100.0;
    double cooling_rate = 0.95;
    int max_iterations = 1000;
    
    // Random number generator for perturbations
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(-5, 5); // Small random perturbation
    
    // Objective function to calculate total overlap error
    auto calculateOverlapError = [&](const vector<Point>& positions) {
        double total_error = 0.0;
        
        // Loop through all adjacent image pairs and calculate overlap error
        for (size_t i = 0; i < positions.size(); ++i) {
            for (size_t j = i + 1; j < positions.size(); ++j) {
                // Calculate overlap error for images i and j
                Rect overlap_region(
                    max(positions[i].x, positions[j].x), 
                    max(positions[i].y, positions[j].y),
                    min(positions[i].x + width, positions[j].x + width) - max(positions[i].x, positions[j].x),
                    min(positions[i].y + height, positions[j].y + height) - max(positions[i].y, positions[j].y)
                );
                
                if (overlap_region.area() > 0) {
                    Mat img1_overlap = images[i](overlap_region);
                    Mat img2_overlap = images[j](overlap_region);
                    total_error += computeMSE(img1_overlap, img2_overlap);
                }
            }
        }
        
        return total_error;
    };

    // Simulated Annealing loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Create a new solution by perturbing current positions
        vector<Point> new_positions = current_positions;
        int idx = rand() % new_positions.size();
        new_positions[idx].x += dis(gen);
        new_positions[idx].y += dis(gen);

        // Calculate the overlap error for the new positions
        double new_error = calculateOverlapError(new_positions);

        // Accept the new positions with a probability based on the temperature
        if (new_error < best_overlap_error || exp((best_overlap_error - new_error) / temperature) > ((double) rand() / RAND_MAX)) {
            current_positions = new_positions;
            if (new_error < best_overlap_error) {
                best_overlap_error = new_error;
                best_positions = new_positions;
            }
        }

        // Cool down the temperature
        temperature *= cooling_rate;
    }

    return best_positions;
}

int main() {
    // Load your images
    vector<Mat> images; // Load 25 images into this vector
    for (int i = 0; i < 25; ++i) {
        images.push_back(imread("image_" + to_string(i) + ".jpg"));
    }

    // Initial positions (rough estimates for the 25 images)
    vector<Point> initial_positions = {{0, 0}, {200, 0}, {400, 0}, {600, 0}, {800, 0}, 
                                       {0, 200}, {200, 200}, {400, 200}, {600, 200}, {800, 200},
                                       {0, 400}, {200, 400}, {400, 400}, {600, 400}, {800, 400},
                                       {0, 600}, {200, 600}, {400, 600}, {600, 600}, {800, 600},
                                       {0, 800}, {200, 800}, {400, 800}, {600, 800}, {800, 800}};
    
    // Assuming all images are of size 200x200
    int width = 200;
    int height = 200;

    // Optimize positions using Simulated Annealing
    vector<Point> optimized_positions = simulatedAnnealing(images, initial_positions, width, height);

    // Print optimized positions
    for (const auto& pos : optimized_positions) {
        cout << "Image position: (" << pos.x << ", " << pos.y << ")" << endl;
    }

    return 0;
}
