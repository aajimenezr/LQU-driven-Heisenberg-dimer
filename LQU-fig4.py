#!/usr/bin/env python
# coding: utf-8

# First order LQU - crytical temperature

# In[4]:


import numpy as np
import matplotlib.pyplot as plt

#parameters

d= .2
alpha= 0.01
J= .5


def Imf(x):
    Imf= 2*d*(J**2 +x**2 + d**2)/(((J + x)**2 + d**2)*((J -x)**2 + d**2))
    return Imf


# ---------------------------------------------------------------------
# Tamaños de fuente (documento a 10pt, RevTeX 4-2 -> caption en \small = 9pt)
# ---------------------------------------------------------------------
FONTSIZE = 10   # etiquetas de ejes (xlabel/ylabel)
TICKSIZE = 9    # números de los ticks

# Ancho de una columna en RevTeX 4-2 (twocolumn, aps, 10pt) ~ 3.4 in
COLWIDTH = 3.4

# NOTA: J, alpha, Imf se asumen ya definidas en tu sesión.

x = np.linspace(0.0001, 5, 10000)
Tc = J / (np.log(np.sqrt((1 + 4 * alpha**2 * (Imf(x))**2) + 2) / np.sqrt(1 + 4 * alpha**2 * (Imf(x))**2)))

fig, ax = plt.subplots(figsize=(COLWIDTH, COLWIDTH * 0.75), dpi=300)
ax.plot(x, Tc)

ax.set_xlim(0, 1.2)
ax.set_xlabel(r"$\omega$", fontsize=FONTSIZE)
ax.set_ylabel(r"Critical temperature $T_c^1$", fontsize=FONTSIZE)
ax.tick_params(axis='both', labelsize=TICKSIZE)
ax.grid(False)   # sin malla, por consistencia con las demás figuras

plt.tight_layout()
plt.savefig('5fig.pdf', dpi=300, bbox_inches='tight')


# In[ ]:




