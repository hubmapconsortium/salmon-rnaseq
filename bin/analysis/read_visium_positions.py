#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from os import walk
import anndata
import manhole
import matplotlib.pyplot as plt
import scanpy as sc

from common import Assay
from plot_utils import new_plot
from typing import Iterable

import sys
import cv2
import tifffile as tf
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from scipy.spatial import distance

TISSUE_COVERAGE_CUTOFF = 0.7

def segment_tissue(img, crop_img, blur_size, morph_kernel_size=153):
    '''
    :param img: RGB Image of the tissue with fiducial spots around it (numpy array 3D)
    :param crop_img:
    :param blur_size:
    :param morph_kernel_size:
    :return:
    '''

    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # gray = fiducial_filter(gray, distance, filter=0)

    # otsu
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    gray = clahe.apply(gray)

    # Apply Gaussian blur to reduce noise
    blurred = cv2.GaussianBlur(gray, (blur_size, blur_size), 0)

    # Apply Otsu's thresholding method to create a binary image
    _, binary = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    # Perform morphological operations to remove small artifacts and fill gaps
    # kernel = np.ones((morph_kernel_size, morph_kernel_size), np.uint8)
    # binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel)
    tissue, _ = fiducial_filter(binary, crop_img, filter=0)
    # connected components to do

    return tissue


def fiducial_filter(img, distance, filter=1):
    '''

    '''

    print(distance)
    print(img.shape)
    distance = (int(distance[0]*img.shape[0]), int(distance[1]*img.shape[1]))
    print(distance)

    if filter:
        mask = np.zeros(img.shape, dtype='uint8')
        mask[:distance[0], :] = 1
        mask[-distance[0]:, :] = 1
        mask[:, :distance[1]] = 1
        mask[:, -distance[1]:] = 1
    else:
        mask = np.ones(img.shape, dtype='uint8')
        mask[:distance[0], :] = 0
        mask[-distance[0]:, :] = 0
        mask[:, :distance[1]] = 0
        mask[:, -distance[1]:] = 0

    return img * mask, distance

def detect_fiducial_spots(img, distance, hough_transform, hough_scale):
    fiducial_crop_img, cropped_pixels = fiducial_filter(img, distance)

    gray = cv2.cvtColor(fiducial_crop_img, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)

    blank_img = np.zeros(blurred.shape)

    # parameters for circle detection
    average_pixels = np.average(cropped_pixels)
    mindist = int(average_pixels / hough_scale[0])
    minr = int(average_pixels / hough_scale[1])
    maxr = int(average_pixels / hough_scale[2])

    # Detect circles using Hough Circle Transform
    circles = cv2.HoughCircles(
        blurred,
        cv2.HOUGH_GRADIENT,
        dp=hough_transform[0],  #best values: [1, 2] - higher there will be more accumulation leading to more circles overlapping but more circles detected
        minDist=mindist, #best values ~ 2x average of min and max radius
        param1=hough_transform[1], # higher of threshold for canny edge detector - greater the value the more strigent the detection of circles
        param2=hough_transform[2], # smaller the more false circles will be detected
        minRadius=minr, #range: 200-300
        maxRadius=maxr #range: 200-300
    )

    if circles is not None:
        circles = np.round(circles[0, :]).astype("int")

        diameters = []
        for (x, y, r) in circles:
            # Draw the circle
            cv2.circle(blank_img, (x, y), r, (255, 0, 0), 2)
            diameters.append(r*2)

        average_diameter = np.average(diameters)

    else:
        print('No Fiducial Spots Found... Exiting')
        sys.exit()

    return circles, average_diameter

    # contour method not as great ##########
    # thresh = cv2.adaptiveThreshold(blurred, 255, cv2.ADAPTIVE_THRESH_MEAN_C,
    #         cv2.THRESH_BINARY_INV, 11, 2)
    #
    # contours, hierarchy = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    # cnt = contours
    #
    # contour_list = []
    # for contour in contours:
    #     approx = cv2.approxPolyDP(contour, 0.01*cv2.arcLength(contour, True), True)
    #     area = cv2.contourArea(contour)
    #     # Filter based on length and area
    #     if (7 < len(approx) < 18) & (900 > area > 200):
    #         # print area
    #         contour_list.append(contour)
    #
    # cv2.drawContours(fiducial_crop_img, contour_list,  -1, (255, 0, 0), 2)
    #######################################


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

    filtered_frames = [frames[i][['X', 'Y']] for i in range(len(frames))]
    normalized_frames = [scaler.fit_transform(frame) for frame in filtered_frames]

    # distance transform - found fiducial spots by frame fiducial spots
    # distance_transform_frames = [distance.cdist(frame, normalized_fiducial_spots) for frame in normalized_frames]
    distance_transform_frames = [distance.cdist(normalized_fiducial_spots, frame) for frame in normalized_frames]
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

        # mask = dist_value[dist_value < threshold]
        # idx = np.where(boolean_mask == 1)
        # dist_sum = dist_value[idx]
        assert (mask.shape[0] == normalized_fiducial_spots.shape[0])

        match_sum.append(sum(mask))
        detected_2_frame.append(thresholded.shape[0] / mask.shape[0])
        frame_2_detected.append(thresholded.shape[0] / distance_transform_frames[0].shape[1])

    # find the frame with the min distance sum and therefore the best frame
    match_idx = np.argmin(match_sum)

    return match_idx, detected_2_frame, frame_2_detected, scaler_fs_detected


def align_N_register(tissue, slide, frame, scaler_fs_detected):
    # normalize min max
    scaler = MinMaxScaler()

    # normalized_tissue = scaler.fit_transform(tissue)
    # resolution = (normalized_tissue.shape[0], normalized_tissue.shape[1])

    inside_diameter = slide['Dia.'].iloc[0]  # should be the same along all columns
    radii = inside_diameter / 2

    # concatenate frame and slide to form block
    filtered_slide = slide[['X', 'Y']]
    filtered_frame = frame[['X', 'Y']]
    filtered_block = pd.concat([filtered_frame, filtered_slide])
    scaler.fit(filtered_block)
    normalized_slide = scaler.transform(filtered_slide)

    # align and scale the normalized slide back to original image resolution
    slide_2_img_res = scaler_fs_detected.inverse_transform(normalized_slide)

    # convert slide of coordinates to image
    new_img = np.zeros(tissue.shape)
    # round to int
    slide_2_img_res = np.round(slide_2_img_res).astype("int")

    # Note: slide_2_img_res idx match exactly that of the index image in new_img
    for i in range(len(slide_2_img_res)):
        # iter_img = np.zeros(tissue.shape)
        # cv2.circle(new_img, (slide_2_img_res[i][0], slide_2_img_res[i][1]), int(radii), (255, 0, 0), -1)
        cv2.circle(new_img, (slide_2_img_res[i][0], slide_2_img_res[i][1]), int(radii), (i + 1, 0, 0), -1)
        #
        # roi = np.where(iter_img > 0)
        # roi_mask = iter_img * tissue
        # tissue_roi = np.where(roi_mask > 0)
        #
        # if np.any(tissue_roi[0]):
        #     fractions.append(tissue_roi[0].shape[0] / roi[0].shape[0])
        # else:
        #     fractions.append(0)

    # convert new_img to int matrix
    new_img = new_img.astype("int")
    # find the percentage of bead covered
    fiducial_idx = np.unique(new_img)[1:]  # do not care about background = 0
    # rois = [np.where(new_img == i) for i in fiducial_idx]
    # rois = [new_img == i for i in fiducial_idx]   #faster but more memory intensive

    fractions = []
    # find fraction of occupancy

    for idx in fiducial_idx:
        roi = new_img == idx
        denom = tissue[roi]  # faster but more memory intensive
        # numerator = tissue[roi[0], roi[1]].flatten()
        numerator = denom[denom > 0].shape[0]
        fractions.append(numerator / denom.shape[0])
    fractions = np.asarray(fractions)

    return fractions

def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath

def get_gpr_df(metadata_dir, img_dir, threshold=None, crop_dim=(0.071313, 0.089899), blur_size=255, morph_kernel_size=153):
    hough_scale = [5.771, 12.545, 11.09]
    hough_parameters = [2, 20, 50]

    gpr_path = list(find_files(metadata_dir, "*.gpr"))[0]
    img_path = list(find_files(img_dir, "*.ome.tiff"))[0]

    gpr = pd.read_table(gpr_path, skiprows=9)

    try:
        img = cv2.imread(str(img_path))
    except:
        img = tf.imread(img_path)

    if threshold is None:
        threshold = gpr['Dia.'].iloc[0] / np.average(img.shape[:2])  # rough estimate

    # big beads = [1, 3, 5, 7] with corresponding inside beads = [2, 4, 6, 8] - # corresponds to block
    # get the frame or big beads that you need for alignment
    frames = []
    tiles = []
    for i in gpr.Block.unique():
        if i % 2:
            # big bead - frame
            frames.append(gpr.loc[gpr['Block'] == i])
        else:
            tiles.append(gpr.loc[gpr['Block'] == i])

    print('Starting tissue segmentation...')
    # segment tissue
    tissue = segment_tissue(img, crop_dim, blur_size=blur_size, morph_kernel_size=morph_kernel_size)
    print('Finish tissue segmentation.')

    print('Starting fiducial spot detection from image...')
    # find fiducial beads in the tissue
    fiducial_spots, fiducial_pixel_diameter = detect_fiducial_spots(img, crop_dim, hough_parameters, hough_scale)
    print('Finish fiducial spot detection.')

    print('Finding best reference frame and slide...')
    # find the right slide
    match_slide_idx, detected_2_frame, frame_2_detected, scaler_fs_detected = slide_match(frames, fiducial_spots,
                                                                                          threshold)
    match_slide = tiles[match_slide_idx].copy()
    match_frame = frames[match_slide_idx].copy()

    fiducial_spatial_diameter = match_frame['Dia.'].iloc[0]
    scale_factor = fiducial_pixel_diameter / fiducial_spatial_diameter

    spot_spatial_diameter = match_slide['Dia.'].iloc[0]

    # align and register tissue with fiducial spots
    fractions = align_N_register(tissue, match_slide, match_frame, scaler_fs_detected)

    # write out a subset of the dataframe for selected capture spots
    # match_slide['Fraction of Tissue Coverage'] = fractions
    match_slide.loc[:, 'Tissue Coverage Fraction'] = fractions
    return match_slide, scale_factor, spot_spatial_diameter

def read_visium_positions(metadata_dir: Path, img_dir: Path, cutoff=0.0):
    gpr_file = list(find_files(metadata_dir, "*.gpr"))[0]
    slide_id = gpr_file.stem
    gpr_df, scale_factor, spot_spatial_diameter = get_gpr_df(metadata_dir, img_dir)
#    gpr_df = pd.read_table(gpr_file, skiprows=9)
#    gpr_df = gpr_df[gpr_df['Block'] == 2]
#    spot_spatial_diameter = gpr_df['Dia.'].iloc[0]
#    scale_factor = 1
    gpr_df = gpr_df.set_index(['Column', 'Row'], inplace=False, drop=True)
    plate_version_number = gpr_file.stem[1]
    barcode_coords_file = Path(f"/opt/data/visium-v{plate_version_number}_coordinates.txt")
    coords_df = pd.read_csv(barcode_coords_file, sep='\t', names=['barcode', 'Column', 'Row'])
    coords_df['Row'] = coords_df['Row'] + 1
    coords_df['Row'] = coords_df['Row'] // 2
    coords_df = coords_df.set_index(['Column', 'Row'])
    gpr_df['barcode'] = coords_df['barcode']
#    gpr_df['Tissue Coverage Fraction'] = 1
    gpr_df = gpr_df[['barcode', 'X', 'Y', 'Tissue Coverage Fraction']]
    gpr_df = gpr_df.reset_index(inplace=False)
    gpr_df = gpr_df.set_index('barcode', inplace=False, drop=True)
    return gpr_df, slide_id, scale_factor, spot_spatial_diameter


def main(metadata_dir: Path, img_dir: Path):
    return read_visium_positions(metadata_dir, img_dir)

if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("metadata_dir", type=Path)
    p.add_argument("img_dir", type=Path)
    args = p.parse_args()

    d = main(args.metadata_dir, args.img_dir)
