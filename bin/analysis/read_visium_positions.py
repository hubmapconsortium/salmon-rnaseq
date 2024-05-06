#!/usr/bin/env python3
import sys
from argparse import ArgumentParser
from os import walk
from pathlib import Path
from typing import Iterable

import cv2
import manhole
import numpy as np
import pandas as pd
import tifffile as tf
from PIL import Image
from scipy.spatial import distance
from scipy.spatial.distance import cdist
from sklearn.preprocessing import MinMaxScaler
import aicsimageio
import ome_utils

def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath

def get_image_from_img_dir(img_dir:Path)->aicsimageio.AICSImage:
    img_files = find_files(img_dir, "*.ome.tif*")
    img_files_list = list(img_files)
    img_file = img_files_list[0]
    img = aicsimageio.AICSImage(img_file)
    return img

def circle_attributes(contour, circularity_threshold=0.85):

    area = cv2.contourArea(contour)
    perimeter = cv2.arcLength(contour, True)

    if perimeter == 0:
        return False, 0, (0, 0)

    circularity = 4 * np.pi * area / (perimeter**2)

    M = cv2.moments(contour)

    if M["m00"] == 0:
        return False, (None,)

    cx = int(M["m10"] / M["m00"])
    cy = int(M["m01"] / M["m00"])

    if circularity > circularity_threshold:
        # Compute the average distance from centroid to contour points
        distances = [cv2.pointPolygonTest(contour, (cx, cy), True)]
        radius = np.mean(distances)

        return True, (cx, cy, abs(radius))

    return False, (None,)


def detect_fiducial_spots_segment_tissue(
    img, threshold, min_neighbors, binary_threshold=125, blur_size=255
):

    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    kernel_close = np.ones((100, 100), np.uint8)  # morphological kernel
    kernel_dilate = np.ones((20, 20), np.uint8)

    # convert image to binary
    _, binary_mask = cv2.threshold(gray, binary_threshold, 255, cv2.THRESH_BINARY)
    binary_mask = cv2.bitwise_not(binary_mask)
    binary_mask = cv2.morphologyEx(binary_mask, cv2.MORPH_CLOSE, kernel_close)
    binary_mask = cv2.dilate(binary_mask, kernel_dilate, iterations=3)

    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edge_img = cv2.Canny(blurred, 50, 200)
    edge_img[binary_mask > 0] = 0

    blank_img = np.zeros(blurred.shape)

    # connected components
    contours, _ = cv2.findContours(edge_img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    new_circles = list(filter(lambda x: x[0], map(circle_attributes, contours)))
    new_circles = np.asarray([x[1] for x in new_circles])
    new_circles[:, 2] = int(np.average(new_circles[:, 2]))
    new_circles = new_circles.astype(int)

    print(f"Number of detected beads BEFORE Outlier Filter: {len(new_circles)}")
    # filter for outliers
    threshold_distance = determine_threshold_distance(new_circles, min_neighbors)
    new_circles = filter_outliers_distance(new_circles, threshold_distance, min_neighbors)
    print(f"Number of detected beads AFTER Outlier Filter: {len(new_circles)}")

    # check for minimum amount of circles detected
    if len(new_circles) < threshold:
        sys.exit(f"Detected fiducial beads: {len(new_circles)} < threshold: {threshold}")

    diameters = []

    for x, y, r in new_circles:
        # Draw the circle
        cv2.circle(blank_img, (x, y), r, (255, 0, 0), 2)
        diameters.append(r * 2)

    average_diameter = np.average(diameters)

    print("Finish fiducial spot detection.")
    print("Starting tissue segmentation...")

    black_pixels = np.where(gray == 0)

    if black_pixels[0].size > 0:
        gray[black_pixels] = 255

    # otsu
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    gray = clahe.apply(gray)

    # Apply Gaussian blur to reduce noise
    blurred = cv2.GaussianBlur(gray, (blur_size, blur_size), 0)

    # Apply Otsu's thresholding method to create a binary image
    _, binary = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    print("Finish tissue segmentation.")

    return new_circles, binary, average_diameter


def slide_match(frames, fiducial_spots, threshold):
    # match percentage for reference spots compared to derived
    match_sum = []

    # min-max normalize fiducial_spots
    scaler_fs_detected = MinMaxScaler()

    scaler_fs_detected.fit(fiducial_spots[:, :2])

    # normalize derived fiducial spots from image
    normalized_fiducial_spots = scaler_fs_detected.transform(fiducial_spots[:, :2])

    # min-max normalize frame fiducial spots
    scaler = MinMaxScaler()
    # filter frames to only important
    # diameter = frames[0]['Dia.'].iloc[0] #should be the same along all columns - 105

    filtered_frames = [frames[i][["X", "Y"]] for i in range(len(frames))]
    normalized_frames = [scaler.fit_transform(frame) for frame in filtered_frames]

    distance_transform_frames = [
        distance.cdist(normalized_fiducial_spots, frame) for frame in normalized_frames
    ]
    match_idx_per_frame = [np.argmin(frame, axis=1) for frame in distance_transform_frames]

    # find percentage of spots found for detected -> frame && frame -> detected fiducial spots
    detected_2_frame = []
    frame_2_detected = []

    for i in range(len(match_idx_per_frame)):
        min_dist_idx = match_idx_per_frame[i]
        dist_value = []
        for j in range(min_dist_idx.shape[0]):
            dist_value.append(distance_transform_frames[i][j, min_dist_idx[j]])

        # threshold
        mask = np.asarray(dist_value)
        thresholded = mask[mask < threshold]

        assert mask.shape[0] == normalized_fiducial_spots.shape[0]

        match_sum.append(sum(mask))
        detected_2_frame.append(thresholded.shape[0] / mask.shape[0])
        frame_2_detected.append(thresholded.shape[0] / distance_transform_frames[0].shape[1])

    # find the frame with the min distance sum and therefore the best frame
    match_idx = np.argmin(match_sum)

    return match_idx, detected_2_frame, frame_2_detected, scaler_fs_detected


def get_rotation_matrix(image):
    # Blur the image slightly to remove noise.
    image = cv2.GaussianBlur(image, (5, 5), 0)
    edges = cv2.Canny(image, 50, 200)

    # Hough Line Transform
    lines = cv2.HoughLines(edges, 1, np.pi / 180, 150)

    if lines is None:
        return None

    angles = []
    for rho, theta in lines[:, 0]:
        # Convert to degrees
        angle = theta * 180 / np.pi - 90
        # Correct angles > 45
        if angle > 45:
            angle -= 180

        # Considering only angles in range [-1, 1]
        if -1 <= angle <= 1:
            angles.append(angle)

    if len(angles) == 0:
        return 0  # return 0 if no small angles detected

    # Return the average of the angles as the rotation
    avg_angles = sum(np.abs(angles)) / len(angles)
    center = (image.shape[1] / 2, image.shape[0] / 2)

    rotation_matrix_2D = cv2.getRotationMatrix2D(center, avg_angles, 1)

    rotation_matrix = np.vstack([rotation_matrix_2D, [0, 0, 1]])

    return rotation_matrix


def filter_outliers_distance(data, threshold_distance, min_neighbors):
    """
    Filters out outliers based on distance.

    - data: Input array of shape (n_samples, n_features).
    - threshold_distance: Maximum distance to consider two points as neighbors.
    - min_neighbors: Minimum number of neighbors a point should have to not be considered an outlier.

    Returns the filtered data.
    """
    num_samples = data.shape[0]
    is_valid = np.zeros(num_samples, dtype=bool)

    for i in range(num_samples):
        distances = np.linalg.norm(
            data - data[i], axis=1
        )  # Compute the distance from the i-th point to every other point.
        num_neighbors = (
            np.sum(distances < threshold_distance) - 1
        )  # Subtract 1 to exclude the point itself.

        # Mark the point as valid if it has the minimum number of neighbors.
        is_valid[i] = num_neighbors >= min_neighbors

    return data[is_valid]


def compute_kth_distances(data, k):
    """
    Compute the distance of each point to its kth nearest neighbor.
    """

    distances = cdist(data, data)
    sorted_distances = np.sort(distances, axis=1)

    # Return the kth distances (k+1 because the 0th distance is the point to itself)
    return sorted_distances[:, k]


def determine_threshold_distance(data, k):
    kth_distances = compute_kth_distances(data, k)
    kth_distances = np.sort(kth_distances)

    # Compute the rate of change between consecutive distances
    derivatives = np.diff(kth_distances)

    # Choose the position where the rate of change is the highest as the "elbow"
    threshold_position = np.argmax(derivatives)

    # Return the kth distance at the "elbow" position
    return kth_distances[threshold_position]


def align_N_register(tissue, slide, frame, scaler_fs_detected, rotational_matrix):
    # normalize min max
    scaler = MinMaxScaler()

    # normalized_tissue = scaler.fit_transform(tissue)
    # resolution = (normalized_tissue.shape[0], normalized_tissue.shape[1])

    inside_diameter = slide["Dia."].iloc[0]  # should be the same along all columns
    radii = inside_diameter / 2

    # concatenate frame and slide to form block
    filtered_slide = slide[["X", "Y"]]
    filtered_frame = frame[["X", "Y"]]
    filtered_block = pd.concat([filtered_frame, filtered_slide])
    scaler.fit(filtered_block)
    normalized_slide = scaler.transform(filtered_slide)

    # align and scale the normalized slide back to original image resolution
    slide_2_img_res = scaler_fs_detected.inverse_transform(normalized_slide)

    scale = scaler.scale_
    min_val = scaler.min_

    # Affine matrix to scale
    affine_matrix_first = [[scale[0], 0, min_val[0]], [0, scale[1], min_val[1]], [0, 0, 1]]

    data_min = scaler_fs_detected.data_min_
    data_max = scaler_fs_detected.data_max_

    # Construct the inverse transformation matrix
    inverse_affine_matrix = [
        [data_max[0] - data_min[0], 0, data_min[0]],
        [0, data_max[1] - data_min[1], data_min[1]],
        [0, 0, 1],
    ]

    affine_matrix = np.dot(inverse_affine_matrix, affine_matrix_first)

    # inverse rotational matrix to transform back to original image space
    r_inv = np.linalg.inv(rotational_matrix)
    # output final transform
    affine_transform = np.dot(r_inv, affine_matrix).T

    # scale back to original resolution
    affine_transform *= 4
    affine_transform[2, 2] = 1
    # convert slide of coordinates to image
    new_img = np.zeros(tissue.shape)
    # round to int
    slide_2_img_res = np.round(slide_2_img_res).astype("int")

    # Note: slide_2_img_res idx match exactly that of the index image in new_img
    for i in range(len(slide_2_img_res)):
        cv2.circle(
            new_img, (slide_2_img_res[i][0], slide_2_img_res[i][1]), int(radii), (i + 1, 0, 0), -1
        )

    # find the percentage of bead covered
    fiducial_idx = np.unique(new_img)[1:]  # do not care about background = 0

    fractions = []
    # find fraction of occupancy
    count = 0

    for idx in fiducial_idx:
        roi = new_img == idx
        denom = tissue[roi]  # faster but more memory intensive
        # numerator = tissue[roi[0], roi[1]].flatten()
        numerator = denom[denom > 0].shape[0]
        fractions.append(numerator / denom.shape[0])
        if numerator / denom.shape[0] > 0.0:
            count += 1

    fractions = np.asarray(fractions)

    return fractions, affine_transform


def downsample_image(image, scale_factor):

    # Resize the image using the new dimensions
    # 5000 by 5000 image was found in testing to provide best tradeoff between runtime and reliable fiducial detection
    resized_image = Image.fromarray(image).resize((5000, 5000), Image.LANCZOS)

    return np.asarray(resized_image)


def get_gpr_df(metadata_dir, img_dir, threshold=None, scale_factor=4, min_neighbors=3):

    gpr_path = list(find_files(metadata_dir, "*.gpr"))[0]
    img_path = list(find_files(img_dir, "*.ome.tiff"))[0]

    gpr = pd.read_table(gpr_path, skiprows=9)

    try:
        img = cv2.imread(str(img_path))
    except:
        img = tf.imread(img_path)
        img = np.transpose(img, (1, 2, 0))

    img = downsample_image(img, scale_factor)

    if threshold is None:
        threshold = gpr["Dia."].iloc[0] / np.average(img.shape[:2])  # rough estimate

    rotational_matrix = get_rotation_matrix(img)
    img = cv2.warpPerspective(img, rotational_matrix, (img.shape[1], img.shape[0]))
    # Does the unrotated image get used for anything after this?

    # big beads = [1, 3, 5, 7] with corresponding inside beads = [2, 4, 6, 8] - # corresponds to block
    # get the frame or big beads that you need for alignment
    frames = []
    tiles = []
    for i in gpr.Block.unique():
        if i % 2:
            # big bead - frame
            frames.append(gpr.loc[gpr["Block"] == i])
        else:
            tiles.append(gpr.loc[gpr["Block"] == i])

    print("Starting fiducial spot detection from image...")
    # find fiducial beads in the tissue
    fiducial_spots, tissue, pixel_diameter = detect_fiducial_spots_segment_tissue(
        img, threshold, min_neighbors
    )
    match_slide_idx, detected_2_frame, frame_2_detected, scaler_fs_detected = slide_match(
        frames, fiducial_spots, threshold
    )
    match_slide = tiles[match_slide_idx].copy()
    match_frame = frames[match_slide_idx].copy()
    fractions, affine_matrix = align_N_register(
        tissue, match_slide, match_frame, scaler_fs_detected, rotational_matrix
    )

    match_slide.loc[:, "Tissue Coverage Fraction"] = fractions

    img = get_image_from_img_dir(img_dir)
    micrometers_per_pixel = ome_utils.get_converted_physical_size(img)['X'].magnitude
    spot_spatial_diameter_micrometers = match_slide["Dia."].iloc[0] #in micrometers
    pixels_per_micrometer = 1 / micrometers_per_pixel #conversion factor from micrometers to pixels
    spot_spatial_diameter_pixels = spot_spatial_diameter_micrometers * pixels_per_micrometer #physical size in pixels

    return match_slide, scale_factor, spot_spatial_diameter_pixels, affine_matrix


def read_visium_positions(metadata_dir: Path, img_dir: Path, cutoff=0.0):
    gpr_file = list(find_files(metadata_dir, "*.gpr"))[0]

    slide_id = gpr_file.stem
    gpr_df, scale_factor, spot_spatial_diameter, affine_matrix = get_gpr_df(metadata_dir, img_dir)

    gpr_df = gpr_df.set_index(["Column", "Row"], inplace=False, drop=True)
    plate_version_number = gpr_file.stem[1]
    barcode_coords_file = Path(f"/opt/data/visium-v{plate_version_number}_coordinates.txt")
    coords_df = pd.read_csv(barcode_coords_file, sep="\t", names=["barcode", "Column", "Row"])
    coords_df["Row"] = coords_df["Row"] + 1
    coords_df["Row"] = coords_df["Row"] // 2
    coords_df = coords_df.set_index(["Column", "Row"])
    gpr_df["barcode"] = coords_df["barcode"]
    gpr_df = gpr_df[["barcode", "X", "Y", "Tissue Coverage Fraction"]]

    gpr_df = gpr_df.reset_index(inplace=False)
    gpr_df = gpr_df.set_index("barcode", inplace=False, drop=True)
    return gpr_df, slide_id, scale_factor, spot_spatial_diameter, affine_matrix


def main(metadata_dir: Path, img_dir: Path):
    return read_visium_positions(metadata_dir, img_dir)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("metadata_dir", type=Path)
    p.add_argument("img_dir", type=Path)
    args = p.parse_args()

    d = main(args.metadata_dir, args.img_dir)
