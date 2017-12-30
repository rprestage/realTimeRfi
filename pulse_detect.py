#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

def pulse_detect(X, Critere, init, return_detections=False):


  N_Ech = X.shape[0];

  ep = 32./2**16;
  #ep = 0.0863;

  Critere_update = 4
  correction     = 1.113110001321999

  N_screen    = 30;
  N_Threshold = 25;

  mu_LP    = np.zeros(N_Ech);
  mu_LP[0] = init;
  LRS3     = np.zeros(3);   # "000" 
  LRS30    = np.zeros(N_screen);  # "000000000000000000000000000000"

  if return_detections:
    WP     = np.zeros(N_Ech);   # Weak pulse
    SP     = np.zeros(N_Ech);   # Strong pulse

  for i in range(1, N_Ech-1):
    # threshold for mean estimate and strong pulse detection
    if X[i]>Critere_update*mu_LP[i-1]:
      LRS3[1:] = LRS3[:-1]
      LRS3[0] = 1
      diff_in = 0;
    else:
      LRS3[1:] = LRS3[:-1]
      LRS3[0] = 0
      diff_in = X[i] - mu_LP[i-1]
    
    # threshold for weak pulse detection
    if X[i]>correction*Critere*mu_LP[i-1]:
      LRS30[1:] = LRS30[:-1]
      LRS30[0] = 1
    else:
      LRS30[1:] = LRS30[:-1]
      LRS30[0] = 0

    
    mu_LP[i] = mu_LP[i-1] + ep * diff_in;
    
    if return_detections:
      if LRS3.sum() == 3:
        SP[i]=1
      if LRS30.sum() > N_Threshold:
        WP[i]=1;


  if return_detections:
    return SP, WP, mu_LP
  else:
    return mu_LP
