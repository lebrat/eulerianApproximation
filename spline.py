import numpy as np
from math import factorial
import matplotlib.pyplot as plt
import os

import sys
from PyQt5.QtWidgets import QApplication, QLabel, QWidget,QLineEdit,QPushButton,QInputDialog,QFormLayout,QMessageBox
from PyQt5.QtGui import QPainter,QPen,QIcon
from PyQt5.QtCore import Qt,QRectF
from PyQt5.QtGui import QMouseEvent,QPainterPath
from scipy.interpolate import spline


class eulerianApproximation():
	def __init__(self,points,m):
		self.controlPoints = points
		# self.points = self.rolledPoints(points,m-1)
		self.points = points
		self.m = m
		self.n = points.shape[1]
		self.nbEul = self.nbEuleriens(self.m,self.n)

	def convol(self,a,b):
		n = a.shape[0]
		result = np.zeros_like(a)
		for i in range(n):
			for j in range(n):
				result[i] += a[j]*b[i-j]
		return result

	def CMCalc(self,mMAX,n):
		T = np.zeros(n)
		T[1] = 1.
		C0 = np.zeros(n)
		C0[0] = 1.
		nbEul = []
		nbEul.append(C0)
		for i in range(1,mMAX):
			CM = np.zeros(n)
			for k,Eul in enumerate(nbEul):			
				CM += 1./float(factorial(i-k))*self.convol(Eul,self.derivExplicit(i-k-1,n))
			nbEul.append(self.convol(CM,T))
			# nbEul.append(CM)

		return nbEul

	def nbEuleriens(self,mMAX,n):
		CM = self.CMCalc(mMAX,n)
		nbEul = []
		for i in range(mMAX):
			nbEul.append(CM[i]*factorial(i))
		return nbEul

	def rolledPoints(self,a,m):
		if m % 2 :
			return .5*(np.roll(a,m//2,axis=1)+np.roll(a,m//2+1,axis=1))
		else:
			return np.roll(a,m//2,axis=1)

	def comb(self,n,k):
		return factorial(n)/(factorial(k)*factorial(n-k))

	def derivExplicit(self,o,n):
		D = np.zeros(n)
		for i in range(o+1):
			D[i] = (-1)**(i)*self.comb(o,i)
		return D
	def computeTheEulerianCurve(self):
		coeffspoints = []
		for k,eul in enumerate(reversed(self.nbEul)):
			Dk = self.derivExplicit(k,self.n)
			coeffspoints += [[self.convol(self.convol(self.points[0,:],Dk),eul)/factorial(self.m-1-k),self.convol(self.convol(self.points[1,:],Dk),eul)/factorial(self.m-1-k)]]
		return coeffspoints


			

class inputNumber(QWidget):
	def __init__(self, parent = None):
		super(inputNumber, self).__init__(parent)
		scriptDir = os.path.dirname(os.path.realpath(__file__))
		self.setWindowIcon(QIcon(scriptDir + os.path.sep + 'eul.png'))
		num,ok = QInputDialog.getInt(self,"integer input dualog","enter a number")
		self.num = int(num)

	
class warningBox(QWidget):
	def __init__(self, parent = None):
		super(warningBox, self).__init__(parent)
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Warning)
		msg.setText("The degree of the appoximation should be lower than the number of control points")
		msg.setWindowTitle("Wrong choice of degree")
		msg.exec_()

class EulerianApproximation(QWidget):
	def __init__(self):
		super().__init__()
		self.initUI()
		# self.setMouseTracking(True)

	def initUI(self):
		self.x = []
		self.y = []
		self.setGeometry(10, 10, 500, 500)
		self.setWindowTitle('EulerianApproximation GUI')
		scriptDir = os.path.dirname(os.path.realpath(__file__))
		self.setWindowIcon(QIcon(scriptDir + os.path.sep + 'eul.png'))
		self.label = QLabel(self)
		# self.label.resize(200, 40)
		self.ready = False
		self.polylineOn = False
		self.interpolationSplineOne = False
		self.eulerianPlot = False
		self.show()



	def mousePressEvent(self, event):
		self.x += [event.x()]
		self.y += [event.y()]
		# print(event.x(),event.y())
		self.ready = True
		self.update()
		
	def keyPressEvent(self, event):
		key = event.key()
		if key == Qt.Key_P :
			print("plotPolyline")
			self.ready = True
			self.polylineOn = True
			self.update()
		if key == Qt.Key_S :
			print("plotSpline")
			self.ready = True
			self.interpolationSplineOne = True
			# computation of the required point
			# level of discretisation
			ld = 200
			nt = np.linspace(0, 1, len(self.x)*ld)
			xtemp = self.x.copy()
			ytemp = self.y.copy()
			xtemp += [xtemp[0]]
			ytemp += [ytemp[0]]
			xtemp += [xtemp[1]]
			ytemp += [ytemp[1]]
			x,y = np.array(xtemp,dtype=float), np.array(ytemp,dtype=float)
			t = np.zeros(x.shape)
			t[1:] = np.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)
			t = np.cumsum(t)
			t /= t[-1]
			timestop = t[-1]
			self.x2 = np.array(spline(t, x, nt)[nt <=timestop])
			self.y2 = np.array(spline(t, y, nt)[nt <=timestop])
	
			siz = self.x2[nt <= t[1]].size	
			time = np.linspace(0,1,siz)
			self.x2[-siz:] = (1-time)*self.x2[-siz:] + time*self.x2[:siz]
			self.y2[-siz:] = (1-time)*self.y2[-siz:] + time*self.y2[:siz]


			self.x2 = np.array(np.round(self.x2),dtype=int)
			self.y2 = np.array(np.round(self.y2),dtype=int)

			self.x2 = self.x2[siz:]
			self.y2 = self.y2[siz:]

			self.update()
		if key == Qt.Key_E :
			self.eulerianPlot = True
			print("plotEulerian")
			numbeASK = inputNumber()
			m = int(numbeASK.num)
			if m > len(self.x):
				war = warningBox()
				raise(ValueError("level of approximation is too large"))
			points = np.array([self.x]+[self.y])
			eA = eulerianApproximation(points,m)
			coeffs = eA.computeTheEulerianCurve()
			nt = 200
			self.x2 = []
			self.y2 = []
			for i in range(len(self.x)):
				t = np.linspace(0,1,nt+1)
				t = t[:-1]
				retx = 0*t
				rety = 0*t
				for deg,coef in enumerate(coeffs):
					retx += coef[0][i]*t**deg/float(factorial(deg))
					rety += coef[1][i]*t**deg/float(factorial(deg))
				self.x2 += [retx]
				self.y2 += [rety]
			self.x2 = np.array(self.x2).ravel()
			self.y2 = np.array(self.y2).ravel()
			self.x2 = np.array(np.round(self.x2),dtype=int)
			self.y2 = np.array(np.round(self.y2),dtype=int)

			self.update()

	def paintEvent(self, e):		
		if self.ready:		
			qp = QPainter()
			qp.begin(self)
			path = QPainterPath()
			for x,y in zip(self.x,self.y):
				path.addRoundedRect(QRectF(x-1, y-1,3,3),1,1)
			qp.fillPath(path,Qt.darkGreen)
			qp.end()
			# print(self.x,self.y)
		if self.polylineOn:
			qp = QPainter()
			qp.begin(self)
			pen = QPen(Qt.darkRed, 2)
			qp.setPen(pen)
			# qp.drawLine(self.x[0],self.y[0],self.x[-1],self.y[-1])
			for i in range(len(self.x)):
				qp.drawLine(self.x[i],self.y[i],self.x[i-1],self.y[i-1])

			path = QPainterPath()
			for x,y in zip(self.x,self.y):
				path.addRoundedRect(QRectF(x-2, y-2,5,5),2,2)
			qp.fillPath(path,Qt.darkGreen)
			qp.end()
		if self.interpolationSplineOne:
			qp = QPainter()
			qp.begin(self)
			pen = QPen(Qt.black, 2)
			qp.setPen(pen)
			for i in range(self.x2.size-1):
				qp.drawLine(self.x2[i],self.y2[i],self.x2[i+1],self.y2[i+1])

			path = QPainterPath()
			for x,y in zip(self.x,self.y):
				path.addRoundedRect(QRectF(x-2, y-2,5,5),2,2)
			qp.fillPath(path,Qt.darkGreen)
			qp.end()

		if self.eulerianPlot:
			qp = QPainter()
			qp.begin(self)
			pen = QPen(Qt.darkCyan, 2)
			qp.setPen(pen)
			for i in range(self.x2.size):
				qp.drawLine(self.x2[i],self.y2[i],self.x2[i-1],self.y2[i-1])

			path = QPainterPath()
			for x,y in zip(self.x,self.y):
				path.addRoundedRect(QRectF(x-2, y-2,5,5),2,2)
			qp.fillPath(path,Qt.darkGreen)
			qp.end()

		self.ready = False
		self.polylineOn = False
		self.interpolationSplineOne = False
		self.eulerianPlot = False







if __name__ == '__main__':
	app = QApplication(sys.argv)
	ex = EulerianApproximation()
	sys.exit(app.exec_())

	# app = QApplication(sys.argv)
	# ex = inputNumber()
	# a = ex.show()
	# print(ex.num)
	# sys.exit(app.exec_())