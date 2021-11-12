# From https://gist.github.com/TengWeiHsu/3d72f1593b03dbcb5bf320b622b3a751
import numpy as np
import math
import time

class OptionPricing:

	def __init__(self,S0,E,T,rf,sigma,iterations):

		self.S0 = S0
		self.E = E
		self.T = T
		self.rf = rf
		self.sigma = sigma
		self.iterations = iterations

	def call_option_simulation(self):

		option_data = np.zeros([self.iterations, 2])
		rand = np.random.normal(0,1, [1, self.iterations])
		stock_price = self.S0*np.exp(self.T*(self.rf - 0.5*self.sigma**2) + self.sigma * np.sqrt(self.T) * rand)
		option_data[:,1] = stock_price - self.E

		# average for the Monte Carlo Method
		average = np.sum(np.amax(option_data, axis = 1))/float(self.iterations)

		return np.exp(-1.0*self.rf*self.T) * average

	def put_option_simulation(self):

		option_data = np.zeros([self.iterations,2])
		rand = np.random.normal(0,1,[1,self.iterations])
		stock_price = self.S0*np.exp(self.T*(self.rf - 0.5*self.sigma**2) + self.sigma * np.sqrt(self.T) * rand)
		option_data[:,1] = self.E - stock_price

		# average for the Monte Carlo Method
		average = np.sum(np.amax(option_data, axis = 1))/float(self.iterations)

		return np.exp(-1.0*self.rf*self.T) * average

if __name__ == "__main__":

	S0=20					#underlying
	E=21				#strike
	T=4/12						#time to maturity
	rf=0.1				#risk-free interest rate
	sigma=0.3			#volatility of the underlying
	iterations=100000000		#number of Monte Carlo iterations

	model = OptionPricing(S0,E,T,rf,sigma,iterations)
	print("Call option price with Monte Carlo approach: ", model.call_option_simulation())
	print("Put option price with Monte Carlo approach: ", model.put_option_simulation())
