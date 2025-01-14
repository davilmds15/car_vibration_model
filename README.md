# car_vibration_model
This code describes the vibrations of a complete car model that passes through an obstacle.

![Repository Logo](image.png)

Car Vibration Model Simulation
This repository contains MATLAB code for simulating the dynamic response of a car model to vibrations. The model considers a combination of translational and rotational dynamics, as well as the effects of suspension and road-induced vibrations.

Overview
The code is designed to simulate the vibrations of a car with the following parameters:

Total car mass and individual wheel masses.
Moments of inertia for yaw and pitch rotations.
Stiffness and damping coefficients for the front and rear suspensions.
Tire stiffness and road profile influences.
Key Features
Simulation of road-induced vibrations on a four-wheel car model.
Time-domain solutions for displacement and rotation of the car body and wheels.
Visualization of dynamic responses for different simulation scenarios.
File Structure
Main MATLAB File: Contains the implementation of the car vibration model, divided into multiple simulations for different road profiles and amplitudes.
Output Plots: Automatically generated to illustrate the dynamic response of the car and wheel displacements during the simulation.
How to Run
Open the .m file in MATLAB.
Define the input parameters for the road profile and car properties at the beginning of the script.
Run the script to generate the plots for the vibration simulations.
Simulation Scenarios
The script includes two main simulations:

Simulation 1: Vibration due to a sinusoidal road profile with a specified amplitude and wavelength for the front axle.
Simulation 2: A different sinusoidal profile affecting only the front-left wheel.
Parameters
You can modify the following parameters in the script:

Mass and Inertia:
m: Total mass of the car.
mf, mr: Masses for front and rear wheels.
Iy, Iz: Moments of inertia.
Stiffness and Damping:
kf, kr: Suspension stiffness (front and rear).
cf, cr: Damping coefficients (front and rear).
ktf, ktr: Tire stiffness.
Road Profile:
A1, A2: Amplitudes of the sinusoidal road profiles.
lambda: Wavelength of the road-induced vibration.
Output
Plots of the dynamic response, including:
Body displacement and rotations (
𝑥
(
𝑡
)
,
𝜙
(
𝑡
)
,
𝜃
(
𝑡
)
x(t),ϕ(t),θ(t)).
Wheel displacements (
𝑥
1
(
𝑡
)
,
𝑥
2
(
𝑡
)
,
𝑥
3
(
𝑡
)
,
𝑥
4
(
𝑡
)
x 
1
​
 (t),x 
2
​
 (t),x 
3
​
 (t),x 
4
​
 (t)).
Road excitation profiles (
𝑦
1
(
𝑡
)
,
𝑦
2
(
𝑡
)
,
𝑦
3
(
𝑡
)
,
𝑦
4
(
𝑡
)
y 
1
​
 (t),y 
2
​
 (t),y 
3
​
 (t),y 
4
​
 (t)).
Author
Name: Davi Lima Mendes dos Santos
University: Universidade Federal de Santa Maria
Course: Dinâmica de Estruturas e Aeroelasticidade
Matricula: 202012282
