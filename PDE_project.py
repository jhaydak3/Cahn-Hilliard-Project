#!/usr/bin/env python 
import numpy as np
import time
import twod_explicit
import twod_ADI


def my_2d_explicit(t0,t_f,x_f,y_f,dx,dy,dt,gamma,D,u0):
	t = np.arange(t0,t_f+dt,step = dt)
	x = np.arange(0,x_f+dx,step = dx)
	y = np.arange(0,y_f+dy,step = dy)
	numx = len(x)
	numy = len(y)
	numt = len(t)
	u0 = np.array(u0,order='F')
	twod_explicit.twod_explicit(u0,numt,dx,dy,dt,gamma,D)
	return x,y,t_f,u0
	
def my_2d_ADI(t0,t_f,x_f,y_f,dx,dy,dt,gamma,D,u0):
	t = np.arange(t0,t_f+dt,step = dt)
	x = np.arange(0,x_f+dx,step = dx)
	y = np.arange(0,y_f+dy,step = dy)
	numx = len(x)
	numy = len(y)
	numt = len(t)
	u0 = np.array(u0,order='F')
	twod_ADI.twod_adi(u0,numt,dx,dy,dt,gamma,D)
	return x,y,t_f,u0
	
def create_initial_rand(numx,numy):
	np.random.seed(seed=34)
	A = -1 + 2*np.random.rand(numx,numy)
	return A
def create_initial_circ(numx,numy):
	A = np.zeros([numx,numy]) - 1
	mdpntx = round(.5*numx)
	mdpnty = round(.5*numy)
	for i in range(0,numx):
		for j in range(0,numy):
			dx = mdpntx - i
			dy = mdpnty - j
			if 2*dx**2 + dy**2 + 5*dx**2*dy - 20*dy**2*dx + 20*dx**3*dy - 10*dy**3*dx < 15:
				A[i,j] = 1
	return A
def main():
	save_str = 'adi_NR_1'
	D = 100
	gamma = .2
	dx = .01; dy = .01
	dt = 1e-11;
	t0 = 0; tf1 = 1e-11
	numSims = 7
	tvec = [0, tf1, tf1*10, tf1*1e2, tf1*1e3, tf1*1e4, tf1*1e5, tf1*1e6]
	time2run = np.zeros(7)
	tout = np.zeros(7)
	x = np.arange(0,1+dx,step = dx)
	numx = len(x)
	#u0 = create_initial_rand(numx,numx)
	u0 = create_initial_circ(numx,numx)
	fh_string = save_str + '_start.txt'
	this_str = 'D = '+ str(D) +', gamma =  ' + str(gamma) + ', dx = ' + str(dx) + ', dy = ' + str(dy)
	this_str = this_str + ', dt = ' + str(dt) + ', t =  ' + str(0) + '\n'
	np.savetxt(fh_string,u0,header = this_str)
	start_time = time.time()
	for i in range(0,numSims):
		#xout,yout,this_t,u0 = my_2d_explicit(tvec[int(i)],tvec[int(i)+1],1,1,dx,dy,dt,gamma,D,u0)
		xout,yout,this_t,u0 = my_2d_ADI(tvec[int(i)],tvec[int(i)+1],1,1,dx,dy,dt,gamma,D,u0)
		time_now = time.time()
		time2run[i] = time_now - start_time
		tvec[int(i)] = this_t
		#write data with file
		fh_string = save_str + '_' + str(i) + '.txt'
		this_str = 'D = '+ str(D) +', gamma =  ' + str(gamma) + ', dx = ' + str(dx) + ', dy = ' + str(dy)
		this_str = this_str + ', dt = ' + str(dt) + ', t =  ' + str(this_t) + ', t_elap = ' + str(time2run[i]) + '\n'
		np.savetxt(fh_string,u0,header = this_str)
		print i
		
if __name__ == "__main__":
	main()