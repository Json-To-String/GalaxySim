import numpy as np
import matplotlib.pyplot as plt


def starAccel(t, vals, mass, xGal, yGal):
	x = vals[1:50] - xGal
	y = vals[51:100] - yGal
	r = np.sqrt(x**2 + y**2)
	a = (-G * mass) / (r**2)
	ax = (a * x) / r
	ay = (a * y) / r

	return([ax, ay])

def galAccel(t, vals, mass, intruder):

	x1 = vals[1:50] - vals[end - 3]
	x2 = vals[1:50] - intruder[end - 3]
	y1 = vals[51:100] - vals[end - 2]
	y2 = vals[51:100] - intruder[end - 2]
	r1 = np.sqrt(x1**2 + y1**2)
	r2 = np.sqrt(x2**2 + y2**2)

	a1 = (-G * mass) / (r1**2)
	a2 = (-G * mass) / (r2**2)

	ax1 = (a1 * x1) / r1
	ax2 = (a2 * x2) / r2
	ay1 = (a1 * y1) / r1
	ay2 = (a2 * y2) / r2

	ax_tot = ax1 + ax2
	ay_tot = ay1 + ay2

	vXGal = vals[end - 1]
	vYGal = vals[end]
	rGal = np.sqrt((vals[end - 3] - intruder[end - 3])**2 + (vals[end - 2] - intruder[end - 2])**2)
	aGal = (-G * mass) / (rGal**2)
	aXGal = aGal * ((vals[end - 3] - intruder[end - 3]) / rGal)
	aYGal = aGal * ((vals[end - 2] - intruder[end - 2]) / rGal)

	return([ax_tot, ay_tot, vXGal, vYGal, aXGal, aYGal])

def galaxies(initCond, mass, direction):

	kpc = 3.0857e19 # kpc to m
	r = np.array([5, 10, 15, 20])
	rAll = kpc * r # m

	theta1 = np.linspace(0, 2*np.pi*r[0] / (r[0]+1), r[0])
	theta2 = np.linspace(0, 2*np.pi*r[1] / (r[1]+1), r[1])
	theta3 = np.linspace(0, 2*np.pi*r[2] / (r[2]+1), r[2])
	theta4 = np.linspace(0, 2*np.pi*r[3] / (r[3]+1), r[3])

	xStars = initCond[0] * kpc + [
				rAll[0] * np.cos(theta1),
				rAll[1] * np.cos(theta2),
				rAll[2] * np.cos(theta3),
				rAll[3] * np.cos(theta4)
			]

	yStars = initCond[1] * kpc + [
				rAll[0] * np.sin(theta1),
				rAll[1] * np.sin(theta2),
				rAll[2] * np.sin(theta3),
				rAll[3] * np.sin(theta4)
			]

	v = lambda r, theta : -theta * np.sqrt((G * mass) / r)

	vX = initCond[3] - direction * [
		v(rAll[1], np.sin(theta1)),
		v(rAll[2], np.sin(theta2)),
		v(rAll[3], np.sin(theta3)),
		v(rAll[4], np.sin(theta4))]

	vY = initCond[4] + direction *[
		v(rAll[1], np.cos(theta1)),
		v(rAll[2], np.cos(theta2)),
		v(rAll[3], np.cos(theta3)),
		v(rAll[4], np.cos(theta4))]

	return([xStars, yStars, vX, vY])

if __name__ == "__main__":

	G =  6.67430e-11
	yr2sec = 31557600 # s
	mSun = 1.989e30 # kg
	mGal = (10e11)*mSun # kg
	kpc = 3.0857e19 # kpc to m

	initCond0 = [0, 0, 0, 0] # kpc kpc m/s m/s
	initCond1 = [25, 175, -100e3, -300e3] # kpc kpc m/s m/s
	initCond2 = [-25, -175, 100e3, 300e3] # kpc kpc m/s m/s

	[xStars0, yStars0, vX0, vY0] = galaxies(initCond0, mGal, -1) # Galaxy on its lonesome

	plt.figure(1)
	plt.plot(xStars0/kpc, yStars0/kpc, 'r*')
	plt.show()
