import math
import random
import matplotlib
import matplotlib.pyplot as plt

atomArray = {}
epsilon = 1.67* 1e-14
epsilon22 = (4/3) * epsilon
epsilon12 = math.sqrt(epsilon * epsilon22)
deltaEm = epsilon12 / 2
sigma = 3.4 * 1e-8
sigma22 = 3.4 * 1e-8
sigma12 = (sigma + sigma22) / 2
omega0 = 1e13
k = 1.38e-16
radius = 3.75e-08
temperature = 1000
nor = 10 # number of rows
noc = nor # number of columns
totalNumOfAtoms = nor * noc
concentration = 25 # %
distance = []
energy = []
plot = 1

set_steps = 0;

class atom:
	neighbor = []
	firstNeighbor = []
	status = 'A'
	locationIndex = 0
	def __init__(self, x, y):
		self.atomX = x
		self.atomY = y
		self.locationIndex = (x - 1) * noc + y

	def getX(self):
		return self.atomX
	def getY(self):
		return self.atomY
	def getIndex(self):
		return self.locationIndex

	def changeStatus(self):
		if self.status == 'A':
			self.status = 'B'
		elif self.status == 'B':
			self.status = 'A'

	def swapeStatus(self, type):
		self.status = type

	def getStatus(self):
		return self.status

	def update(self, atomArray):
		self.neighbor = []
		self.firstNeighbor = []
		indexList = []
		x = self.atomX
		y = self.atomY

		if x == 1:
			if y != 1 and y != nor:
				indexList = [(x-1)*nor+y+1,(x)*nor+y,(x)*nor+y+1,(x-1)*nor+y-1,(x)*nor+y-1,nor*(nor-1)+y-1,nor*(nor-1)+y,nor*(nor-1)+y+1]
			elif y == 1:
				indexList = [(x-1)*nor+y+1,(x)*nor+y,(x)*nor+y+1,(x)*nor+y-1,(x+1)*nor+y-1,nor*(nor-1)+y,nor*(nor-1)+y+1,nor*nor+y-1]
			elif y == nor:
				indexList = [1,(x-1)*nor+y-1,(x-1)*nor+y+1,(x)*nor+y-1,(x)*nor+y,nor*(nor-1)+1,nor*(nor-1)+y-1,nor*(nor-1)+y]
		elif x == nor:
			if y != 1 and y != nor:
				indexList = [y-1,y,y+1,(x-2)*nor+y-1,(x-2)*nor+y,(x-2)*nor+y+1,(x-1)*nor+y-1,(x-1)*nor+y+1]
			elif y == 1:
				indexList = [y,y+1,nor,(x-2)*nor+y,(x-2)*nor+y+1,(x-1)*nor+y-1,(x-1)*nor+y+1,nor*nor+y-1]
			elif y == nor:
				indexList = [1,y-1,y,(x-3)*nor+y+1,(x-2)*nor+y-1,(x-2)*nor+y,(x-2)*nor+y+1,(x-1)*nor+y-1]
		elif y == 1:
			indexList = [(x-2)*nor+y,(x-2)*nor+y+1,(x-1)*nor+y-1,(x-1)*nor+y+1,(x)*nor+y-1,(x)*nor+y,(x)*nor+y+1,(x+1)*nor+y-1]
		elif y == nor:
			indexList = [(x-3)*nor+y+1,(x-2)*nor+y-1,(x-2)*nor+y,(x-2)*nor+y+1,(x-1)*nor+y-1,(x-1)*nor+y+1,(x)*nor+y-1,(x)*nor+y]
		else:
			indexList = [(x-2)*nor+y-1,(x-2)*nor+y,(x-2)*nor+y+1,(x-1)*nor+y-1,(x-1)*nor+y+1,(x)*nor+y-1,(x)*nor+y,(x)*nor+y+1]

		for i in indexList:
			self.neighbor.append(atomArray[i])

		L = nor
		for atom in self.neighbor:
			spanX = abs(atom.getX() - self.atomX)
			spanY = abs(atom.getY() - self.atomY)
			if (spanX > L / 2):
				spanX = L - abs(atom.getX() - self.atomX)
			if (spanY > L / 2):
				spanY = L - abs(atom.getY() - self.atomY)
			if ((spanX + spanY) <= 1):
				if atom.getX() != self.atomX or atom.getY() != self.atomY:
					self.firstNeighbor.append(atom)

	def update_N2(self, atomArray):
		self.neighbor = []
		L = nor
		for i in range(1,totalNumOfAtoms+1):
			spanX = abs(atomArray[i].getX() - self.atomX)
			spanY = abs(atomArray[i].getY() - self.atomY)
			if (spanX > L / 2):
				spanX = L - abs(atomArray[i].getX() - self.atomX)
			if (spanY > L / 2):
				spanY = L - abs(atomArray[i].getY() - self.atomY)
			if (spanX <= 1) and (spanY <= 1):
				if atomArray[i].getX() != self.atomX or atomArray[i].getY() != self.atomY:
					self.neighbor.append(atomArray[i])

	def getNeighbor(self):
		return self.neighbor

	def get1stNeighbor(self):
		return self.firstNeighbor


	def getEnergy(self, distance):
		self.energy = 0
		L = nor

		for atoms in self.neighbor:
			spanX = abs(atoms.getX() - self.atomX)
			spanY = abs(atoms.getY() - self.atomY)
			if (spanX > L / 2):
				spanX = L - abs(atoms.getX() - self.atomX)
			if (spanY > L / 2):
				spanY = L - abs(atoms.getY() - self.atomY)

			self.atomDist = distance * (spanX **2 + spanY **2) ** (1/2.0)

			if self.status == 'A' and atoms.getStatus() == 'A':
				self.energy += 4 * epsilon * ((sigma / (self.atomDist)) ** 12 - ((sigma / (self.atomDist)) ** 6))
			elif self.status == 'B' and atoms.getStatus() == 'B':
				self.energy += 4 * epsilon22 * ((sigma22 / (self.atomDist)) ** 12 - ((sigma22 / (self.atomDist)) ** 6))
			else:
				self.energy += 4 * epsilon12 * ((sigma12 / (self.atomDist)) ** 12 - ((sigma12 / (self.atomDist)) ** 6))

		return self.energy / 2

	def getEnergy_Simple(self, distance):
		self.energy = 0

		self.firstD = 4 * epsilon * ((sigma / (self.atomDist)) ** 12 - ((sigma / (self.atomDist)) ** 6))
		self.SecondD = 4 * epsilon * ((sigma / (self.atomDist)) ** 12 - ((sigma / (self.atomDist)) ** 6))

		self.energy = self.firstD + self.SecondD
		return self.energy

	# override toString
	def __str__(self):
		return ("(" + self.status + ", " + str(self.atomX) + ", " + str(self.atomY) + ")")
	def __unicode__(self):
		return u("(" + self.status + ", " + str(self.atomX) + ", " + str(self.atomY) + ")")
	def __repr__(self):
		return ("(" + self.status + ", " + str(self.atomX) + ", " + str(self.atomY) + ")")

def energy_v_distance():
	for dist in range(300,600):
		calDist = dist * 1e-10
		distance.append(calDist) # x axies

		energyAtDist = 0
		for i in range (1,totalNumOfAtoms+1):
			energyAtDist += atomArray[i].getEnergy(calDist)

		energy.append(energyAtDist)

	print(min(energy))
	print(energy.index(min(energy)))
	print(distance[energy.index(min(energy))])	

	# plot
	plt.plot(distance, energy)
	plt.xlabel('Distance')
	plt.ylabel('Energy')
	plt.show()

def MMC():
	energy = []
	count = 0
	accept = 0
	numOfRun = []

	totalEnergy = 0

	for i in range(1,totalNumOfAtoms+1):
		totalEnergy += (atomArray[i].getEnergy(radius))

	while count < 50000:
		count += 1
		# energy of state i
		randNum1 = random.randint(0,numOfBatom-1) # random B atom
		randNum2 = random.randint(1,totalNumOfAtoms)      # random atom
		#print(randNum1, randNum2)
		
		energyAtStatei = atomArray[BatomList[randNum1]].getEnergy(radius) + atomArray[randNum2].getEnergy(radius)

		# swap
		if atomArray[randNum2].getStatus() != 'B':
			atomArray[BatomList[randNum1]].changeStatus()
			atomArray[randNum2].changeStatus()
		elif atomArray[randNum2].getStatus() == 'B':
			continue

		# energy of state j
		energyAtStatej = atomArray[BatomList[randNum1]].getEnergy(radius) + atomArray[randNum2].getEnergy(radius)

		deltaH = energyAtStatej - energyAtStatei
		#print('delta H:', deltaH)
		if deltaH >= 0:
			xi = random.random()
			if temperature == 0:
				atomArray[BatomList[randNum1]].changeStatus()
				atomArray[randNum2].changeStatus()

			elif xi >= math.exp(-deltaH / (k * temperature)):
				atomArray[BatomList[randNum1]].changeStatus()
				atomArray[randNum2].changeStatus()

			elif xi < math.exp(-deltaH / (k * temperature)):
				if randNum2 not in BatomList:
					BatomList[randNum1] = randNum2
				numOfRun.append(accept+1)
				totalEnergy += deltaH
				energy.append(totalEnergy)

				accept += 1
				#print(accept)
		elif deltaH < 0:
			if randNum2 not in BatomList:
				BatomList[randNum1] = randNum2
			numOfRun.append(accept+1)
			totalEnergy += deltaH
			energy.append(totalEnergy)

			accept += 1
			#print(accept)
	print('Total accept: ', accept)

	plt.plot(numOfRun, energy)
	plt.xlabel('Steps')
	plt.ylabel('Energy')

def KMC():
	t = 0
	Time = []
	energy = []
	deltaH_list = []
	totalEnergy = 0
	count = 0

	for i in range(1,totalNumOfAtoms+1):
		totalEnergy += (atomArray[i].getEnergy(radius))

	while count < 10000:
		Rn = 0
		deltaH_list = []
		total_event = []
		count += 1
		for Batom in BatomList:
			for atom in atomArray[Batom].get1stNeighbor():
				energyAtStatei = atomArray[Batom].getEnergy(radius) + atom.getEnergy(radius)

				# swape
				keep = atomArray[Batom].getStatus()
				atomArray[Batom].swapeStatus(atom.getStatus())
				atom.swapeStatus(keep)

				energyAtStatej = atomArray[Batom].getEnergy(radius) + atom.getEnergy(radius)

				# swape back
				keep = atomArray[Batom].getStatus()
				atomArray[Batom].swapeStatus(atom.getStatus())
				atom.swapeStatus(keep)

				deltaH = energyAtStatej - energyAtStatei
				deltaE = deltaH + deltaEm

				deltaH_list.append(deltaH)

				if temperature == 0:
					omega_event = 0
				else:
					omega_event = omega0 * math.exp(deltaE / (- k * temperature))

				total_event.append(omega_event)

		for i in total_event:
			Rn += i

		if Rn != 0:
			total_event = [i / Rn for i in total_event]

		xi1 = 1 - random.random()
		cul = 0
		target_index = 0
		for ind, value in enumerate(total_event):
			cul += value
			if cul >= xi1 and (cul - value) < xi1:
				target_index = ind
				break

		atomInd = target_index // 4
		direction = target_index % 4
		#　call the "direction atom" inside the neighbot list of the target atom
		atom_direction = atomArray[BatomList[atomInd]].get1stNeighbor()[direction]
		# update B atom list
		if atom_direction.getStatus() == 'A':
			BatomList[atomInd] = atom_direction.getIndex()
			
		# swap
		keep = atomArray[BatomList[atomInd]].getStatus()
		atomArray[BatomList[atomInd]].swapeStatus(atom_direction.getStatus())
		atom_direction.swapeStatus(keep)

		# Energy change
		totalEnergy += deltaH_list[target_index]
		energy.append(totalEnergy)

		# delta_T
		xi2 = 1 - random.random()

		if Rn != 0:
			delta_T = math.log(1 / xi2) / Rn
		else:
			delta_T = float('inf')

		t += delta_T
		Time.append(t)
		
		plt.subplot(1,3,3)
		plt.plot(Time, energy)
		plt.xlabel('Time')
		plt.ylabel('Energy')
		

def plot_atoms():
	plt.subplot(1, 3, plot)
	for i in range (1,totalNumOfAtoms+1):
		if atomArray[i].getStatus() == 'A':
			plt.plot(atomArray[i].getX(),atomArray[i].getY(),'bo', markersize = 20)
		if atomArray[i].getStatus() == 'B':
			plt.plot(atomArray[i].getX(),atomArray[i].getY(),'ro', markersize = 20)

def double_check():
	forA = 0
	forB = 0
	for i in range (1,totalNumOfAtoms+1):
		if atomArray[i].getStatus() == 'A':
			forA += 1
		elif atomArray[i].getStatus() == 'B':
			forB += 1
	print("A: ", forA, "B: ", forB)
	print("B list length:", len(BatomList))

# Build atoms matrix
for i in range (1,nor+1):
	for j in range (1,noc+1):
		newatom = atom(i,j)
		atomArray.setdefault((i - 1) * noc + j, newatom)
# Build neighbor list
for i in range (1,totalNumOfAtoms+1):
	atomArray[i].update(atomArray)

reportNum = 1
neighbor = atomArray[reportNum].getNeighbor()
print('neighbor list of atom {0}: {1}'.format(atomArray[reportNum].getIndex(), neighbor))

first_neighbor = atomArray[reportNum].get1stNeighbor()
print ('first neighbor list:',first_neighbor)

# Randomize B atom position
numOfBatom = int(totalNumOfAtoms * concentration / 100)
BatomList = random.sample(range(1,totalNumOfAtoms+1),numOfBatom)
#print(BatomList)

for i in BatomList:
	atomArray[i].changeStatus()
	#print(atomArray[i])

# lowest energy is -8.224e-10 at 3.75e-08
# energy_v_distance()


if concentration == 0:
	numOfBatom = 1  # for index exception of 0 concentration 

double_check()

plot_atoms() # original ditribution
plt.title('Original Distribution')
plot += 1
# MMC()  # Metropolis Monte Carlo
KMC()  # Kinetic Monte Carlo


# varificaiton
double_check()
plot_atoms()  # plot atoms' distribution
plt.title('KMC')
'''
for i in BatomList:
	print(atomArray[i])
'''
plt.show()