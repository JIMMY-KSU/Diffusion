import numpy as np

class Diffusion2D():
		
	def __init__(self,t,Nx,Ny):
		self.a = 0; self.b = 1
		self.c = 0; self.d = 1
		self.D = 1 #Diffusion constant
		self.t = t; self.Nx = Nx; self.Ny = Ny;
		
		self.Nt = 0.1	#Time iterations

		self.hx = (self.b - self.a)/float((self.Nx -1))	#meshsize
		self.hy = (self.b - self.a)/float((self.Ny -1))	#meshsize
		self.ht = t/float(self.Nt-1)

	def MC2D_U(self,walkers):
		
		l0 =  sqrt(2*self.D*self.ht)
		x = np.zeros(walkers)
		y = np.zeros(walkers)

		for k in range(100):
			dx = random.randint(1,5,len(x))
			dy = random.rand(len(y))
			#b = nonzero(dx<0.25)
			#c = nonzero(dy<0.25)

			r_dx = array( dx==1, dtype=int)
			x[:] = x[:] + r_dx[:]*l0

			r_dy = array( dx ==2, dtype=int)
			y[:] = y[:] + r_dy[:]*l0
		
			l_dx = array (dx ==3,dtype=int)
			x[:] = x[:] - l_dx[:]*l0
			
			l_dy = array (dx== 4,dtype=int)
			y[:] = y[:] - l_dy[:]*l0
			
			#y = y.tolist()
			#x = x.tolist()	
		
			for i in range(len(x)):
				if x[i] <0:
					x[i] =0
				if x[i] >1:
					x[i] =0
			for j in range(len(y)):
				if y[i] <0:
					y[i] = 0
				if y[i] >1:
					y[i] = 0

		plot(x,y,'o')
		axis([0,1,0,1])

	def explicit(self):

		Nt,Nx,Ny,hx,hy = self.Nt,self.Nx, self.Ny, self.hx, self.hy

		if self.ht > (1./2)*(((hx*hy)**2)/(hx**2 + hy**2)):
			print " The numerical scheme is not stable for this condition"

		u = np.zeros([Nx,Ny])
		unew = np.zeros([Nx,Ny])

		#Boundary conditions
		u[0,:] 	= 0
		u[-1,:] = 0
		u[:,0] 	= 1
		u[:,-1] = 0
		unew[:,:] = u[:,:]
	
		for t in range(1,Nt+1):
			unew[1:-1,1:-1] = u[1:-1,1:-1] + self.ht*((u[2:,1:-1] - 4*u[1:-1,1:-1] + u[:-2,1:-1] +u[1:-1,2:] + \
			u[1:-1,:-2])/self.hy**2)

			u[:,:] = unew[:,:]
		return u

	def Jacobi(self,max_iterations):
		
		Nt,Nx,Ny = self.Nt,self.Nx, self.Ny
		alpha = float(self.ht)/(self.hx*self.hx)
		u = np.zeros([Nx,Ny])
		ug = np.zeros([Nx,Ny])
		unew = np.zeros([Nx,Ny])

		#Boundary conditions
		u[0,:] 	= 0
		u[-1,:] = 0
		u[:,0] 	= 1
		u[:,-1] = 0
		unew[:,:] = u[:,:]
		ug[:,:] = u[:,:]

		iterations = 0
		diff = 1

		#Vectorized code
		for t in range(1,self.Nt):
			u[:,:] = unew[:,:]
			iterations = 0
			diff = 1

			while iterations < max_iterations and diff > 1E-6:
				diff = 0
				unew[1:-1,1:-1] = (1/(1+4*alpha))*(alpha*(ug[2:,1:-1] + ug[:-2,1:-1] + ug[1:-1,2:] + ug[1:-1,:-2]) + u[1:-1,1:-1])
				diff = sum(unew[:,:]-ug[:,:])
				ug[:,:] = unew[:,:]
				iterations += 1

		return unew

	def exact(self):

		x = np.linspace(self.a, self.b,self.Nx)
		y = np.linspace(self.c, self.d, self.Ny)
		u_exact = np.zeros([self.Nx,self.Ny])

		for m in range(1,(self.Nx*2)):
			for n in range(1,(self.Nx*2)):
				u_exact = u_exact + self.fmn(x,y,m,n,self.t)
		
		return u_exact

	def fmn(self,x,y,m,n,t):

		Amn = (-np.cos(m*np.pi) + 1)/(m*np.pi)*(-2/(n*np.pi))
		return np.sin(n*np.pi*x)*np.sin(m*np.pi*y)*1#*np.exp(np.pi**2*(m**2 + n**2)*t)*Amn

	def plotting(self,u):

		x = np.linspace(0,1,self.Nx)
		y = np.linspace(0,1,self.Ny)
		x ,y = np.meshgrid(x,y)

		import matplotlib.pylab as plt
		from matplotlib import cm
		from mpl_toolkits.mplot3d import Axes3D

		fig = plt.figure()
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(x, y, u, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
		fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()

def execute():
	d1 = Diffusion2D(t=0.001,Nx=100,Ny=100)
	u = d1.exact()
	#u = d1.explicit()
	#u = d1.Jacobi(500)
	#u = d1.MC2D_U(1000)
	d1.plotting(u)

if __name__ == '__main__':
	execute()
