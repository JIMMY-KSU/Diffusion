from scitools.std import *

class Diffusion:
	"""
	This class solves the diffusion equation in 1D by means of the following algorithms: explicit,implicit,
	cranck-Nicholson -scheme, and by consider a random walker in 1D.
	"""
	def __init__(self,t):
		self.a = 0
		self.b = 1
		self.Nx = 36
		self.hx = (self.b - self.a)/float((self.Nx -1))
		self.Nt = 50 	
		self.t = t
		self.ht = t/float(self.Nt-1)
		self.D = 0.1
		self.alpha = (self.D*self.ht)/(self.hx*self.hx)

	def exact(self):
		N = 100
		x = linspace(self.a, self.b,self.Nx)
		u_exact = zeros(self.Nx)
		tic = time.time()
		for n in range(1,N+1):
			An = -2/(n*pi)
			u_exact[:] = u_exact[:] + An*sin(n*pi*x[:])*exp(-self.D*((n*pi)**2)*self.t)

		u_exact[:] = 1 - x[:] + u_exact[:]
		t_exact = linspace(self.a,self.b, len(x))
		toc = time.time()
		return t_exact, u_exact

	def explicit(self):
		Nx,Nt,t, alpha = self.Nx, self.Nt, self.t, self.alpha
		if alpha >= 0.5:
			print "The explicit scheme is unstable for this spatio-temporal gridding"
			sys.exit(1)

		print alpha
		u = zeros(Nx)
		unew = zeros(Nx)
		u[0] = 1 ; unew[0] = 1; u[-1] = 0; unew[-1] = 0
		j = 1
		while j <= (Nt-1):
			unew[1:-1] = alpha*u[:-2] + (1-2*alpha)*u[1:-1] + alpha*u[2:]
			u[:] = unew[:]
			j +=1

		return unew
	
	def implicit(self):

		Nx, alpha, Nt = self.Nx, self.alpha, self.Nt
		b = (1+2*alpha); a = (-1*alpha); c = a
		unew = zeros(Nx)
		u = zeros(Nx)
		vnew = zeros(Nx)
		v = zeros(Nx)
		grid = linspace(self.a, self.b,Nx)

		#initial condition
		u[0] = 1 ; unew[0] = 1; u[-1] = 0; unew[-1] = 0
		g = lambda x:1-x
		v = u - g(grid)
		vnew = unew - g(grid)
		
		for i in range(1,Nt):
			new = self.tridiagonal(a,b,c,vnew,v,Nx)
			v[:] = vnew[:]
	
		unew = vnew + g(grid)
		return unew

	def cranck_nicolson(self):

		Nx, alpha, Nt = self.Nx, self.alpha, self.Nt

		b = (2+2*alpha); a = (-1*alpha); c = a		
		unew = zeros(Nx)
		u = zeros(Nx)
		temp = zeros(Nx)
		v = zeros(Nx)
		vnew = zeros(Nx)
		grid = linspace(self.a, self.b,Nx)
		
		#Initial condition
		u[0] = 1 ; unew[0] = 1; u[-1] = 0; unew[-1] = 0
		g = lambda x:1-x
		v = u - g(grid)
		vnew = unew - g(grid)
		
		#Time loop
		for i in range(1,Nt):
			for j in range(1,len(v)-1):
				temp[j] = alpha*v[j-1] + alpha*v[j+1] + (2-2*alpha)*v[j]

			vnew = self.tridiagonal(a,b,c,vnew,v,Nx)
			v[:] =vnew[:]

		unew = vnew + g(grid)
		return unew

	def tridiagonal(self, a,b,c,vnew,v,Nx):
			
		temp = zeros(Nx)
		btemp = b
		vnew[1] = v[1]/btemp

		#forward substitution
		for i in range(1,self.Nx-1):
			temp[i] = c/btemp
			btemp = b-a*temp[i]
			vnew[i] = (v[i] - a*vnew[i-1])/btemp

		#backward substitution
		for i in range(self.Nx-2,0,-1):
			vnew[i] = vnew[i] - temp[i+1]*vnew[i+1]
		
		return vnew

	def MC_uniform(self,walkers):

		l0 =  sqrt(2*self.D*self.ht)
		temp_1 = int((1./l0)) + 1
		position = zeros(temp_1,int)
		temp = zeros(temp_1,int)
		position[0] = walkers
		temp[0] = walkers

		for j in range(self.Nt):
			for i in range(len(position)):
				if temp[i] != 0:
					for k in range(temp[i]):
						r = random.random()
						#right
						if r > 0.5 and i!= 0 and i!= (len(position)-2):
							position[i+1] += 1
							position[i] -= 1
						if r > 0.5 and i== 0:
							position[i+1] += 1
						if r > 0.5 and i == (len(position)-2):
							position[i+1] = 0

						#left
						if r < 0.5 and i!= 0 and i!= 1:
							position[i-1] +=1
							position[i] -=1
							
						if r< 0.5 and i==1:
							position[i] -=1	
			temp[:]= position[:]

		position_scale = (1./position[0])*position
		t_uniform = linspace(self.a, self.b, len(position))
		return t_uniform, position_scale

	def MC_normal(self,walkers):
		Nbins = 10
		t = linspace(0,1,Nbins)
		dudes = []
		
		for i in range(walkers):
			dudes.append(0)
		total = walkers
		count = 0

		for j in range(self.Nt):
			count = 0
			total = len(dudes)
			while count < total:
				epsilon = random.normal(0,(1/sqrt(2)))
				dudes[count] = dudes[count] + sqrt(2*self.D*self.ht)*epsilon
				if dudes[count]< 0:
					del dudes[count]
					dudes.append(0)
				if dudes[count] >1:
					del dudes[count]
		
				count += 1	

		
		hist, bin_edges = histogram(dudes)
		hist = (1./max(hist)) * hist
		return hist, t
		

def error():
	Nx = 10
	explicit_error = zeros(Nx)
	implicit_error = zeros(Nx)
	MC_uniform_error = zeros(Nx)
	for i in range(1,11):
		d1 = Diffusion(t=0.3)
		t_exact,u_exact = d1.exact()
		u_explicit = d1.explicit()	
		u_implicit = d1.implicit()

def compare_MC_exact():
	max_error = []
	sum_error = []
	walkers = []
	for i in range(500,20000,500):
		d1 = Diffusion(t=0.2)
		t_exact,u_exact = d1.exact()
		t_uniform , mc_uniform = d1.MC_uniform(i)
		print len(u_exact)
		print len(mc_uniform)
		diff = u_exact[:] - mc_uniform[:]
		temp = max(abs(diff))
		max_error.append(temp)
		temp = sum(diff)
		sum_error.append(temp)
		temp = i
		print temp
		walkers.append(temp)

	from matplotlib.pylab import plt
	plt.figure(2)	
	plt.plot(walkers,max_error,'o-')
	plt.xlabel('Number of walkers')
	plt.ylabel('Maximum error')
	plt.savefig('mcuniform_error1.eps')

	plt.figure(3)
	plt.plot(walkers,sum_error,'o-')
	plt.xlabel('Number of walkers')
	plt.ylabel('Accumulated error')
	plt.savefig('mcuniform_error2.eps')

def execute():
	d1 = Diffusion(t=1)
	#t_exact,u_exact = d1.exact()
	hist,Nbins = d1.MC_normal(100)
	#plot(t_exact,u_exact,'*-')
	#hold('on')
	plot(hist)

if __name__ == '__main__':
	execute()
	#compare_MC_exact()
