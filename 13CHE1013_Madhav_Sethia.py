{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    
    "## Fitting Kinetics\n",
    "\n",
    "## Introduction\n",
    "We have been given an assignment to improve the productivity of a particular process.\n",
    "$$\n",
    "    A + B \\rightarrow C   \n",
    "$$\n",
    "$$\n",
    "    A + C \\rightarrow D\n",
    "$$\n",
    "$$\n",
    "   B \\rightarrow degradation\n",
    "$$\n",
    "The first two reactions are irreversible and of 1st order wrt each reactant.  The last reaction is <i>second</i> order wrt B.  The kinetic constants are $k_1, k_2, k_3$."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "If the initial quantity of A is $N_{A0}$, write down differential equations with boundary conditions that govern the conversion of A, B, C and D.\n",
    "\n",
    "$$\n",
    "    N_A(t+\\Delta t) - N_A(t) = V(-r_1 - r_2)\\Delta t \\rightarrow \\frac{dN_A}{dt} = V(-r_1-r_2)\n",
    "$$\n",
    "Similarly:\n",
    "$$\n",
    "    \\frac{dN_B}{dt} = V(-r_1-r_3) + R\n",
    "$$\n",
    "$$\n",
    "    \\frac{dN_C}{dt} = V(r_1 - r_2)\n",
    "$$\n",
    "$$\n",
    "    \\frac{dN_D}{dt} = V(r_2)\n",
    "$$\n",
    "$$\n",
    "    \\frac{dV}{dt} = \\frac{R}{\\rho}\n",
    "$$\n",
    "The boundary conditions are: \n",
    "$$\n",
    "    \\left[ N_A, N_B, N_C, N_D \\right]_{t=0} = \\left[ N_{A0}, 0, 0, 0 \\right]\n",
    "$$\n",
    "Where, $r_1$, $r_2$, $r_3$ are the volumetric rates of the three reactions.\n",
    "$$\n",
    "    r_1 = k_1C_AC_B\n",
    "$$\n",
    "$$\n",
    "    r_2 = k_2C_AC_C\n",
    "$$\n",
    "$$\n",
    "    r_3 = k_3C_B^2\n",
    "$$\n",
    "And, $\\left[ C_A, C_B, C_C \\right] = \\left[ \\frac{N_A}{V}, \\frac{N_B}{V}, \\frac{N_C}{V}  \\right]$."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "\n",
    "Using data from the file \"ExamProblemData.csv\" fit the kinetic constants $k_1, k_2, k_3$ for each temperature.  The data headers are as follows:  Col1, Col2, Col3 are the time (s), $C_A, C_D$ for T=250K.  Col4, Col5, Col6 are for T=300K.  Col7, Col8, Col9 are for T=350K and Col10, Col11, Col12 are for T=400K.  Concentrations are in kmol/m3, $N_{A0} = 100 kmol$ and $R = 1$ kmol/$m^3$. \n",
    "\n",
    "\n",
    
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 253,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "import scipy\n",
    "import pandas as pd\n",
    "import seaborn as sn\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Hence Col1 is to be renamed t_250, Col2 will be C_A_250, Col3 C_D_250 etc.  "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 254,
   "metadata": {
    "collapsed": false
   },
   "outputs": [],
   "source": [
    "headers = ['t_250', 'C_A_250', 'C_D_250', 't_300', 'C_A_300', 'C_D_300', 't_350', 'C_A_350', 'C_D_350', 't_400', 'C_A_400', 'C_D_400']\n",
    "df = pd.read_csv(\"ExamProblemData.csv\",  #The name of the csv file\n",
    "                 header = 0,            
    "                 names = headers         
    "                 )     "
   ]
  },
  