import numpy as np
import cv2
import random
import math
import os

# Function to compute Mean Squared Error between two images
def compute_mse(img1, img2):
    diff = cv2.absdiff(img1, img2)
    diff = diff.astype(np.float32)
    diff = diff ** 2
    mse = np.sum(diff) / (img1.size)
    return mse

# Simulated Annealing for optimizing the positions
def simulated_annealing(images, initial_positions, width, height):
    current_positions = initial_positions.copy()
    best_positions = initial_positions.copy()
    best_overlap_error = float('inf')
    
    # Parameters for simulated annealing
    temperature = 100.0
    cooling_rate = 0.95
    max_iterations = 1000
    
    # Objective function to calculate total overlap error
    def calculate_overlap_error(positions):
        total_error = 0.0
        
        # Loop through all adjacent image pairs and calculate overlap error
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                # Calculate overlap region for images i and j
                overlap_x1 = max(positions[i][0], positions[j][0])
                overlap_y1 = max(positions[i][1], positions[j][1])
                overlap_x2 = min(positions[i][0] + width, positions[j][0] + width)
                overlap_y2 = min(positions[i][1] + height, positions[j][1] + height)
                
                overlap_width = overlap_x2 - overlap_x1
                overlap_height = overlap_y2 - overlap_y1
                
                if overlap_width > 0 and overlap_height > 0:
                    # Extract the overlapping region from both images
                    img1_overlap = images[i][overlap_y1:overlap_y2, overlap_x1:overlap_x2]
                    img2_overlap = images[j][overlap_y1:overlap_y2, overlap_x1:overlap_x2]
                    
                    # Compute MSE between the overlapping regions
                    total_error += compute_mse(img1_overlap, img2_overlap)
        
        return total_error

    # Simulated Annealing loop
    for iter in range(max_iterations):
        # Create a new solution by perturbing current positions
        new_positions = current_positions.copy()
        idx = random.randint(0, len(new_positions) - 1)
        new_positions[idx] = (new_positions[idx][0] + random.randint(-5, 5), new_positions[idx][1] + random.randint(-5, 5))
        
        # Calculate the overlap error for the new positions
        new_error = calculate_overlap_error(new_positions)
        
        # Accept the new positions with a probability based on the temperature
        if new_error < best_overlap_error or math.exp((best_overlap_error - new_error) / temperature) > random.random():
            current_positions = new_positions
            if new_error < best_overlap_error:
                best_overlap_error = new_error
                best_positions = new_positions
        
        # Cool down the temperature
        temperature *= cooling_rate
    
    return best_positions

# Main function to execute the tiling and optimization
def main():
    # Load your images (for example, 25 images)
    images = []
    for i in range(25):
        img = cv2.imread(f'image_{i}.jpg')
        images.append(img)
    
    # Initial positions (rough estimates for the 25 images)
    initial_positions = [(0, 0), (200, 0), (400, 0), (600, 0), (800, 0), 
                         (0, 200), (200, 200), (400, 200), (600, 200), (800, 200),
                         (0, 400), (200, 400), (400, 400), (600, 400), (800, 400),
                         (0, 600), (200, 600), (400, 600), (600, 600), (800, 600),
                         (0, 800), (200, 800), (400, 800), (600, 800), (800, 800)]
    
    # Assuming all images are of size 200x200
    width = 200
    height = 200
    
    # Optimize positions using Simulated Annealing
    optimized_positions = simulated_annealing(images, initial_positions, width, height)
    
    # Print optimized positions
    for pos in optimized_positions:
        print(f'Optimized image position: {pos}')

if __name__ == '__main__':
    main()
