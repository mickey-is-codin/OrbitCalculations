#https://github.com/smit2300/OrbitCalculations

import math
import random
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d


class GetInput(object):

	def get_eye(self):
		self.eye_radius = float(input("Radius of the eyeball in mm: "))		#Method for getting eyeball radius input
		return self.eye_radius

	def get_bleb(self):		
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

	#****CHANGE INDEXING TO START AT LAST POINT****#


	desired_volumes = [10, 50, 100, 150, 200, 250, 300]
	eye_radius = 11.55
	bleb_radius = [0 for i in range(0,len(desired_volumes))]
	tolerance = 0.005
	dx = 0.0001
	rev_volume = [0, 10, 50, 100, 150, 200, 250]
	bleb_volume = [0 for i in range(0,len(desired_volumes))]
	bleb_area = [0 for i in range(0,len(desired_volumes))]
	meet_point_x = [0 for i in range(0,len(desired_volumes))]
	meet_point_y = [0 for i in range(0,len(desired_volumes))]
	quad1_area = [0 for i in range(0,len(desired_volumes))]
	rev_volume = [0 for i in range(0,len(desired_volumes))]
	surface_instance = [0 for i in range(0,len(desired_volumes))]
	surface_area_bleb = [0 for i in range(0,len(desired_volumes))]
	surface_area_eye = [0 for i in range(0,len(desired_volumes))]

	print("\n*--------------------------------*\n")

	for i in range(0,len(desired_volumes)):
		while rev_volume[i] < (desired_volumes[i] - tolerance):
			eye_volume = get_volume(eye_radius)																		#Eyeball volume in mm^3
			bleb_volume[i] = get_volume(bleb_radius[i])																#Bleb volume in mm^3	

			eye_area = get_area(eye_radius)																			#Eyeball surface area in mm^2
			bleb_area[i] = get_area(bleb_radius[i])																	#Bleb surface area in mm^2
			
			#Bleb Equation: y = sqrt(r^2 - x^2) or x = sqrt(r^2 - y^2)
			#Eyeball Equation: y = R - sqrt(R^2 - x^2) or x = sqrt((-1)*(y-2*R)y)

			meet_point_x[i] = intersection_x(bleb_radius[i], eye_radius)											#The x location of the intersection of the two circles
			meet_point_y[i] = intersection_y(bleb_radius[i], eye_radius)											#The y locaiton of the intersection of the two circles
			
			quad1_area[i] = quad1_x_section(bleb_radius[i], eye_radius, meet_point_x[i], dx)						#The cross sectional area between the two circles contained in quadrant I
			
			rev_volume[i] = rev_integral(bleb_radius[i], eye_radius, meet_point_x[i], meet_point_y[i], dx)			#Volume of the cross section revolved around the y axis
			
			surface_instance[i] = SurfaceIntegral(bleb_radius[i], eye_radius, meet_point_x[i], meet_point_y[i], dx)	
			surface_area_bleb[i] = surface_instance[i].surface_bleb()
			surface_area_eye[i] = surface_instance[i].surface_eye()	

			bleb_radius[i] = bleb_radius[i] + .01

			if i < len(desired_volumes)-1:
				bleb_volume[i+1] = bleb_volume[i]
				bleb_area[i+1] = bleb_area[i]
				meet_point_x[i+1] = meet_point_x[i]
				meet_point_y[i+1] = meet_point_y[i]
				quad1_area[i+1] = quad1_area[i]
				rev_volume[i+1] = rev_volume[i]
				surface_instance[i+1] = surface_instance[i]
				surface_area_bleb[i+1] = surface_instance[i]
				surface_area_eye[i+1] = surface_area_eye[i]
				bleb_radius[i+1] = bleb_radius[i]

		print("Intersection Point of Bleb and Eye: ", meet_point_x[i], meet_point_y[i])
		print("Bleb Radius: %.2f" % bleb_radius[i])
		print("Volume of Bleb: %.2f" % rev_volume[i])	
		print("Surface Area of Bleb: %.2f" % surface_area_bleb[i])
		print("Surface Area of Eye: %.2f" % surface_area_eye[i])
		print("\n*--------------------------------*\n")

	PlotData(bleb_radius, rev_volume, surface_area_bleb, surface_area_eye)
	Draw2DData(bleb_radius, eye_radius, rev_volume, surface_area_bleb, surface_area_eye)
	print("\n*--------------------------------*\n")
	Statistics(bleb_radius, eye_radius, rev_volume, surface_area_bleb, surface_area_eye)
	#Draw3DData(bleb_radius, eye_radius, rev_volume, surface_area_bleb, surface_area_eye)

def get_volume(r):
	volume = (4/3) * math.pi * r**3																					#Using the generic volume of a sphere equation
	return volume

def get_area(r):
	area = 4 * math.pi * r**2																						#Using the generic surface area of a sphere equation
	return area

def intersection_x(r, R):	
	meet_point_x = math.sqrt((-1)*r**4 + 4 * r**2 * R**2)/(2 * R)													#Solve both circle equations for y and set them equal. Wolfram Alpha used for algebra.
	return meet_point_x

def intersection_y(r, R):
	meet_point_y = (r**2) / (2 * R)																					#Solve both circle equations for x and set them equal. Wolfram Alpha used for algebra.
	return meet_point_y

def quad1_x_section(r, R, meet_point, dx):

	N = int(meet_point / dx)																						#Number of points along x axis that we will use as solutions of integration
	nums = [x*dx for x in range(0,N)]																				#Setting an array to iterate through for all of the solution points of integration

	eye_points = []																									#Blank array for all of the points solved along the eyeball's curve
	bleb_points = []																								#Blank array for all of the points solved along the bleb's curve
	area = 0																										#Will build the area by iteratively multiplying the point by dx to give incremental solution

	for n in range(0,N-1):
		eye_points.append(R - math.sqrt(R**2 - nums[n]**2))
		bleb_points.append(math.sqrt(r**2 - nums[n]**2))
		area = area + bleb_points[n]*dx - eye_points[n]*dx

	return area

def rev_integral(r, R, meet_point_x, meet_point_y, dx):

	N_1 = int(meet_point_y / dx)																					#Using two different ranges because of two integration bounds
	N_2 = int((r - meet_point_y) / dx)																				#First range is 0 to meet_y, second is meet_y to r because we integrate wrt y

	nums_1 = [y*dx for y in range(0,N_1)]																			#Create an array of all the y-axis solution points
	nums_2 = [y*dx + meet_point_y for y in range(0,N_2)]															#Create another array of the y-axis solution points for the eyeball curve

	solutions = []
	eye_solutions = []																								#Blank array for all of the points solved along the eyeball's curve
	bleb_solutions = []																								#Blank array for all of the points solved along the bleb's curve
	volume_eye = 0																									#Will build the integral by iteratively multiplying the point by dx to give incremental solution
	volume_bleb = 0																									#Will build the integral by iteratively multiplying the point by dx to give incremental solution
	rev_volume = 0
	volume = 0

	for n in range(0,len(nums_1)):
		eye_solutions.append(math.pi * (nums_1[n]*(2*R-nums_1[n])))
		volume_eye = volume_eye + eye_solutions[n]*dx

	for n in range(0,len(nums_2)):
		bleb_solutions.append(math.pi * (r**2 - nums_2[n]**2))
		volume_bleb = volume_bleb + bleb_solutions[n]*dx
	
	rev_volume = volume_eye + volume_bleb
	return rev_volume


class SurfaceIntegral(object):

	def __init__(self, r, R, meet_point_x, meet_point_y, dx):
		self.r = r
		self.R = R
		self.meet_point_x = meet_point_x
		self.meet_point_y = meet_point_y
		self.dx = dx

	def surface_bleb(self):

		self.N_2 = int((self.r - self.meet_point_y) / self.dx)
		self.nums_2 = [y*self.dx + self.meet_point_y for y in range(0,self.N_2)]

		self.bleb_solutions = []
		self.surface_area_bleb = 0

		for n in range(1,len(self.nums_2)):
			self.f_prime_y = (-1) * (self.nums_2[n]) / math.sqrt((self.r**2 - self.nums_2[n]**2))
			self.ds = math.sqrt(1 + self.f_prime_y**2)
			self.f_y = math.sqrt(self.r**2 - self.nums_2[n]**2)

			self.bleb_solutions.append(2 * math.pi * self.f_y * self.ds)	
			self.surface_area_bleb = self.surface_area_bleb + self.bleb_solutions[n-1]*self.dx

		return self.surface_area_bleb

	def surface_eye(self):

		self.N_1 = int(self.meet_point_y / self.dx)
		self.nums_1 = [y*self.dx for y in range(0,self.N_1)]

		self.eye_solutions = []
		self.surface_area_eye = 0

		for n in range(1,len(self.nums_1)):
			self.f_prime_y = (self.R - self.nums_1[n]) / (math.sqrt(self.nums_1[n]) * math.sqrt(2*self.R - self.nums_1[n]))
			self.ds = math.sqrt(1 + self.f_prime_y**2)
			self.f_y = math.sqrt(self.nums_1[n]) * math.sqrt(2*self.R - self.nums_1[n])

			self.eye_solutions.append(2 * math.pi * self.f_y * self.ds)
			self.surface_area_eye = self.surface_area_eye + self.eye_solutions[n-1]*self.dx

		return self.surface_area_eye

def PlotData(r, vol, area_bleb, area_eye):

	fig1 = plt.figure(figsize=(10,10))
	ax = fig1.add_subplot(1, 1, 1)

	plt.title("Surface Area of Bleb Compared with Surface Area of Eye for Given Volumes", fontsize=16)
	plt.xlabel("Volume of Bleb (mm^3)")
	plt.ylabel("Surface Area (mm^2)")

	plt.scatter(
		x=vol,
		y=area_bleb,
		marker='o',
		c='r',
		label='Surfacea Area of Bleb'
	)

	plt.scatter(
		x=vol,
		y=area_eye,
		marker='x',
		c='b',
		label='Outer Coverage Surface Area'
	)

	plt.legend()

	plt.show()
	fig1.savefig('surface_area_comparison.png', dpi=100)

def Draw2DData(r, R, vol, area_bleb, area_eye):

	#colors = ['b', 'k', 'r', 'cyan', 'purple', 'green', 'w']
	r_color = lambda: random.uniform(0.0,1.0)
	colors = [(r_color(),r_color(),r_color()) for i in range(0,len(r))]
	blebs = [0 for x in range(0,len(colors))]

	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(1, 1, 1)

	eyeball = plt.Circle(
		(0,11.55),
		11.55,
		color='k',
		fill=False
	)

	sclera = plt.Circle(
		(0,11.55),
		2.5,
		color='g',
		fill=True,
		alpha=0.5
	)

	pupil = plt.Circle(
		(0,11.55),
		1.0,
		color='k',
		fill=True
	)

	for i in range(0,len(colors)):
		blebs[i] = plt.Circle(
			(0,0),
			r[i],
			fill=False,
			color=colors[i],
			alpha=0.3
		)

	ax.set_xlim((-15, 15))
	ax.set_ylim((0, 30))

	ax.add_artist(eyeball)
	ax.add_artist(sclera)
	ax.add_artist(pupil)

	for i in range(0,len(blebs)):
		ax.add_artist(blebs[i])
		#ax.fill_between(blebs[i], eyeball, where=blebs[i]>eyeball, facecolor=colors[i], alpha=0.1)

	plt.show()

def Statistics(r, R, vol, area_bleb, area_eye):
	#print("Bleb Radii: ", r)
	#print("Bleb Volume: ", vol)
	#print("Bleb Surface Area: ", area_bleb)
	#print("Eye Surface Area: ", area_eye)
	print("Statistical Analysis")

	radius_volume_ratio = []
	radius_bleb_area_ratio = []
	radius_eye_area_ratio = []

	for i in range(0,len(r)):
		radius_volume_ratio.append(r[i] / vol[i])
		radius_bleb_area_ratio.append(r[i] / area_bleb[i])
		radius_eye_area_ratio.append(r[i] / area_eye[i])

	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(1, 1, 1)

	plt.title("Surface Area of Bleb Compared with Surface Area of Eye for Given Volumes", fontsize=16)
	plt.xlabel("Volume of Bleb (mm^3)")
	plt.ylabel("Surface Area (mm^2)")

	plt.scatter(
		x=vol,
		y=radius_volume_ratio,
		marker='o',
		c='r',
		label='Bleb Radius : Bleb Volume'
	)
	plt.plot(
		vol,
		radius_volume_ratio,
		'-'
		'r',
	)

	plt.scatter(
		x=vol,
		y=radius_bleb_area_ratio,
		marker='x',
		c='b',
		label='Bleb Radius : Bleb Area'
	)
	plt.plot(
		vol,
		radius_bleb_area_ratio,
		'-'
		'b',
	)

	plt.scatter(
		x=vol,
		y=radius_eye_area_ratio,
		marker='^',
		c='g',
		label='Bleb Radius : Eye Area'
	)
	plt.plot(
		vol,
		radius_eye_area_ratio,
		'-'
		'g',
	)

	plt.legend()

	plt.show()


def Draw3DData(r, R, vol, area_bleb, area_eye):

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	r_color = lambda: random.uniform(0.0,1.0)
	colors = [(r_color(),r_color(),r_color()) for i in range(0,len(r))]

	u = np.linspace(0, 2 * np.pi, 100)
	v = np.linspace(0, np.pi, 100)
	x = 11.55 * np.outer(np.cos(u), np.sin(v))
	y = 11.55 * np.outer(np.sin(u), np.sin(v)) + 11.55
	z = 11.55 * np.outer(np.ones(np.size(u)), np.cos(v))

	ax.plot_surface(x, y, z, color='gray', alpha=0.2)

	for i in range(0,len(r)):
		u = np.linspace(0, 2 * np.pi, 100)
		v = np.linspace(0, np.pi, 100)
		x = r[i] * np.outer(np.cos(u), np.sin(v))
		y = r[i] * np.outer(np.sin(u), np.sin(v))
		z = r[i] * np.outer(np.ones(np.size(u)), np.cos(v))

		ax.plot_surface(x, y, z, color=colors[i], alpha=0.2)

	plt.show()

if __name__ == '__main__':
	main()
