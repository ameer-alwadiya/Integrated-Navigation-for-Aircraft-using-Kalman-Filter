# Integrated Navigation System with Extended Kalman Filter

This project involves the design and implementation of an integrated navigation system that combines GPS, IMU, and air-data inputs. The system utilizes the Extended Kalman Filter (EKF) to estimate 12 states, including position, velocity, attitude, and wind components. The design process includes deriving Jacobian matrices, implementing the EKF, removing sensor biases, detecting and diagnosing faults, and addressing cyber attack events.

## Project Structure

- `SID230255137.m`: MATLAB script file for implementing the Extended Kalman Filter and related tasks.
- `dataTask4.mat`: Data file used for testing and validation.

## Overview

### 1. Kalman Filter Design

**1.1. Jacobian Matrices**
Derivation of Jacobian matrices necessary for the EKF with a nonlinear system.

- **F Matrix Elements**:
  - `F(5, 4) = -r_m`
  - `F(6, 7) = -(981*cos(theta)*sin(phi))/100`

- **G Matrix Elements**:
  - `G(4, 6) = -v`
  - `G(9, 5) = -sin(phi)/cos(theta)`

- **H Matrix Elements**:
  - `H(7, 7) = 1`
  - `H(11, 4) = -w/(u^2*(w^2/u^2 + 1))`

**1.2. Implementing EKF**
- Estimations of 12 states and trajectories are plotted.
- Observations include initial spikes and convergence to stable trajectories.

### 2. Removing Sensor Biases

**2.1. Estimation Comparison**
- Comparison between biased and unbiased measurements.
- Observed differences and the effects of biases on state estimations.

**2.2. Removing Biases**
- Incorporation of bias terms as variables.
- Improved estimation of states after including biases.

### 3. Sensor Fault Detection and Diagnosis

**3.1. Identifying Sensor Faults**
- Detection of faults through measurement innovations.
- Analysis of fluctuations indicating fault occurrences.

**3.2. Removing Sensor Faults**
- Inclusion of bias terms as random walk processes.
- Improved state estimations and fault detection.

### 4. Dealing with Cyber Attack Events

**4.1. Cyber Attack Detection**
- Detection of attacks on angle of attack sensor.
- Use of innovation and threshold adjustments for detecting anomalies.

**4.2. Mitigating Cyber Attack Effects**
- Adjustment of standard deviation to handle increased uncertainty.
- Comparison of estimations before and after adjustments.

## How to Run

1. **Setup**: Ensure you have MATLAB installed with access to the required toolboxes.
2. **Load Data**: Load the `dataTask4.mat` file into MATLAB.
3. **Run Script**: Execute the `SID230255137.m` script to perform the Kalman Filter design, estimation, and analysis tasks.
4. **View Results**: Check the generated figures for insights into state estimations, bias removal, fault detection, and cyber attack mitigation.

## Requirements

- MATLAB R2021a or later
- Signal Processing Toolbox
- Control System Toolbox

## Contact

For any questions or issues, please contact ameer.alwadiya@outlook.com.
