#https://github.com/smit2300/OrbitCalculations

import math
import numpy as np

class GetInput(object):

	def get_eye(self):
		self.eye_radius = float(input("Radius of the eyeball in mm: "))		#Method for getting eyeball radius input
		return self.eye_radius

	def get_bleb(self):
		#self.bleb_radius = float(input("Radius of the bleb in mm: "))		#Method for getting bleb radius input
		
		self.bleb_radius = input("List of radii of the bleb in mm separated by spaces: ").split(' ')

		for i in range(0,len(self.bleb_radius)):
			self.bleb_radius[i] = float(self.bleb_radius[i])
		
		return self.bleb_radius

	def get_dx(self):
		self.dx = float(input("Space division for integration (dx), recommended 0.001 or smaller: "))		#Method for getting space division input
		return self.dx

def main():
	#----Start Main Function----#

	#----General Notation: R = Eyeball radius, r = bleb_radius----#

	fetch_input = GetInput()									#Create an instance of the input getting class so we can use it for further input

	desired_volumes = [10, 50, 100, 150, 200, 250, 300]
	eye_radius = 11.55
	bleb_radius = [0 for i in range(0,len(desired_volumes))]
	tolerance = 0.005
	dx = 0.0001
	rev_volume = [0 for i in range(0,len(desired_volumes))]
	bleb_volume = [0 for i in range(0,len(desired_volumes))]
	bleb_area = [0 for i in range(0,len(desired_volumes))]
	meet_point_x = [0 for i in range(0,len(desired_volumes))]
	meet_point_y = [0 for i in range(0,len(desired_volumes))]
	quad1_area = [0 for i in range(0,len(desired_volumes))]
	rev_volume = [0 for i in range(0,len(desired_volumes))]
	surface_area = [0 for i in range(0,len(desired_volumes))]

	for i in range(0,len(desired_volumes)):
		while rev_volume[i] < (desired_volumes[i] - tolerance):
			eye_volume = get_volume(eye_radius)							#Eyeball volume in mm^3
			bleb_volume[i] = get_volume(bleb_radius[i])						#Bleb volume in mm^3

			eye_area = get_area(eye_radius)								#Eyeball surface area in mm^2
			bleb_area[i] = get_area(bleb_radius[i])							#Bleb surface area in mm^2

			#Bleb Equation: y = sqrt(r^2 - x^2) or x = sqrt(r^2 - y^2)
			#Eyeball Equation: x = sqrt(y) * sqrt(2*R - y)

			meet_point_x[i] = intersection_x(bleb_radius[i], eye_radius)			#The x location of the intersection of the two circles
			meet_point_y[i] = intersection_y(bleb_radius[i], eye_radius)			#The y locaiton of the intersection of the two circles

			quad1_area[i] = quad1_x_section(bleb_radius[i], eye_radius, meet_point_x[i], dx)			#The cross sectional area between the two circles contained in quadrant I

			rev_volume[i] = rev_integral(bleb_radius[i], eye_radius, meet_point_x[i], meet_point_y[i], dx)	#Volume of the cross section revolved around the y axis
			surface_area[i] = surface_integral(bleb_radius[i], eye_radius, meet_point_x[i], meet_point_y[i], dx)	#Surface of the volume created by revolution

			bleb_radius[i] = bleb_radius[i] + .01

		print("Bleb Radius: ", bleb_radius[i])
		print("Volume: ", rev_volume[i])

def get_volume(r):
	volume = (4/3) * math.pi * r**3									#Using the generic volume of a sphere equation
	return volume

def get_area(r):
	area = 4 * math.pi * r**2										#Using the generic surface area of a sphere equation
	return area

def intersection_x(r, R):	
	meet_point = math.sqrt((-1)*r**4 + 4 * r**2 * R**2)/(2 * R)		#Solve both circle equations for y and set them equal. Wolfram Alpha used for algebra.
	return meet_point

def intersection_y(r, R):
	meet_point = r**2 / (2 * R)										#Solve both circle equations for x and set them equal. Wolfram Alpha used for algebra.
	return meet_point

def quad1_x_section(r, R, meet_point, dx):

	N = int(meet_point / dx)										#Number of points along x axis that we will use as solutions of integration
	nums = [x*dx for x in range(0,N)]								#Setting an array to iterate through for all of the solution points of integration

	eye_points = []													#Blank array for all of the points solved along the eyeball's curve
	bleb_points = []												#Blank array for all of the points solved along the bleb's curve
	area = 0														#Will build the area by iteratively multiplying the point by dx to give incremental solution

	for n in range(0,N-1):
		eye_points.append(R - math.sqrt(R**2 - nums[n]**2))
		bleb_points.append(math.sqrt(r**2 - nums[n]**2))
		area = area + bleb_points[n]*dx - eye_points[n]*dx

	#print("Cross sectional area in the first quadrant: %.2f square mm" % area)
	#print("\n-----------------------------------------------------------------\n")

	return area

def rev_integral(r, R, meet_point_x, meet_point_y, dx):
	
	#----Volume of the shape made by both curves rotated about y-axis----#

	'''Made with the understanding that 
		V = integral from 0 to meet_point_y of pi * bleb equation^2 
		+ integral from meet_point_y to r of pi * eye equation^2'''

	N_1 = int(meet_point_y / dx)									#Using two different ranges because of two integration bounds
	N_2 = int((r - meet_point_y) / dx)								#First range is 0 to meet_y, second is meet_y to r because we integrate wrt y

	nums_1 = [y*dx for y in range(0,N_1)]							#Create an array of all the y-axis solution points
	nums_2 = [y*dx for y in range(0,N_2)]							#Create another array of the y-axis solution points for the eyeball curve

	eye_integral = []												#Blank array for all of the points solved along the eyeball's curve
	bleb_integral = []												#Blank array for all of the points solved along the bleb's curve
	volume_eye = 0													#Will build the integral by iteratively multiplying the point by dx to give incremental solution
	volume_bleb = 0													#Will build the integral by iteratively multiplying the point by dx to give incremental solution

	for n in range(0,N_2-1):
		eye_integral.append(math.pi * (math.sqrt(nums_2[n]) * math.sqrt(2 * R - nums_2[n]))**2)
		volume_eye = volume_eye + eye_integral[n]*dx

	for n in range(0,N_1-1):
		bleb_integral.append(math.pi * (math.sqrt(r**2 - nums_1[n]**2))**2)
		volume_bleb = volume_bleb + bleb_integral[n]*dx

	rev_volume = volume_eye + volume_bleb							#Sum of the two integrals shoudl be the volume
	#print("Volume of the bleb contained in the eye (displacement of the eye): %.2f cubic mm" % rev_volume)
	#print("\n-----------------------------------------------------------------\n")

	return rev_volume

def surface_integral(r, R, meet_point_x, meet_point_y, dx):

	#----Surface area of the shape made by both curves rotated about y-axis----#

	'''Made with the understanding that 
		Area = integral from 0 to meet_point_y of 2 * pi * bleb equation
		+ integral from meet_point_y to r of 2 * pi * eye equation'''

	#----Function is the same as that used to gather volume of revolution values, but using 2 * pi * radius to drive the equations----#


	N_1 = int(meet_point_y / dx)
	N_2 = int((r - meet_point_y) / dx)

	nums_1 = [y*dx for y in range(0,N_1)]
	nums_2 = [y*dx for y in range(0,N_2)]

	eye_integral = []
	bleb_integral = []
	surface_eye = 0
	surface_bleb = 0

	for n in range(0,N_2-1):
		eye_integral.append(2 * math.pi * (math.sqrt(nums_2[n]) * math.sqrt(2 * R - nums_2[n])))
		surface_eye = surface_eye + eye_integral[n]*dx


	for n in range(0,N_1-1):
		bleb_integral.append(2 * math.pi * (math.sqrt(r**2 - nums_1[n]**2)))
		surface_bleb = surface_bleb + bleb_integral[n]*dx

	surface_area = surface_eye + surface_bleb
	#print("Full surface area of the bleb contained in the eye: %.2f square mm" % surface_area)
	#print("Surface area of the bleb receiving drug: %.2f square mm" % surface_eye)
	#print("Surface area of the eye that the bleb is underneath: %.2f square mm" % surface_bleb)
	#print("\n-----------------------------------------------------------------\n")

	return surface_area


if __name__ == '__main__':
	main()
