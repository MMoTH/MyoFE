# -*- coding: utf-8 -*-
"""
Created on Thu July 31 18:58:51 2022

@author: Hossein
"""

import numpy as np
import json

from scipy.integrate import odeint

class growth():

    def __init__(self, growth_structure,
                parent_circulation):
