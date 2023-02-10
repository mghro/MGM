#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/10/23 3:39 PM

@author: alejandrobertolet
"""
import numpy as np
from scipy.stats import poisson, gamma

class MicrodosimetricGammaModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetFunctionsAndParameters(pars)

    def getBDD(self, yF):
        return self.BDD(yF, *self.BDD_pars)

    def getBDI(self, yF):
        return self.BDI(yF, *self.BDI_pars)

    def getSBD(self, yF):
        return self.SBD(yF, *self.SBD_pars)

    def getSBI(self, yF):
        return self.SBI(yF, *self.SBI_pars)

    def getN_sites(self, yF):
        return self.N_sites(yF, *self.N_sites_pars)

    def getN_sites_with_DSB(self, yF):
        return self.N_sites_with_DSB(yF, *self.N_sites_with_DSB_pars)

    def getGamma_par1(self, yF):
        return self.gamma_par1(yF, *self.gamma_par1_pars)

    def getGamma_par2(self, yF):
        return self.gamma_par2(yF, *self.gamma_par2_pars)

    def getGamma(self, yF):
        return self.gamma_func(yF, self.getGamma_par1(yF), self.getGamma_par2(yF))

    def SetFunctionsAndParameters(self, pars=None):
        if pars is None:
            self.BDD = self.linear_func
            self.BDD_pars = np.array([1.15261696])
            self.BDI = self.exp_func
            self.BDI_pars = np.array([9.24573603e+02, 4.20536461e-03])
            self.SBD = self.linear_func
            self.SBD_pars = np.array([0.96784241])
            self.SBI = self.exp_func
            self.SBI_pars = np.array([1.52868524e+02, 8.72750460e-03])
            self.N_sites = self.linear_plus_exp_func
            self.N_sites_pars = np.array([4.20443273e-01, 5.67814904e+02, 1.04384649e-02])
            self.N_sites_with_DSB = self.linquad_func
            self.N_sites_with_DSB_pars = np.array([0.13781355, 0.00060011])
            self.gamma_par1 = self.quadratic_func
            self.gamma_par1_pars = np.array([8.65049809e-05, 5.65761513e-03, 1.35721323e+00])
            self.gamma_par2 = self.quadratic_func
            self.gamma_par2_pars = np.array([-6.62740411e-05,  1.15279062e-03,  1.55850260e+00])
        else:
            if 'BDD' in pars:
                self.BDD = pars['BDD'][0]
                self.BDD_pars = pars['BDD'][1]
            if 'BDI' in pars:
                self.BDI = pars['BDI'][0]
                self.BDI_pars = pars['BDI'][1]
            if 'SBD' in pars:
                self.SBD = pars['SBD'][0]
                self.SBD_pars = pars['SBD'][1]
            if 'SBI' in pars:
                self.SBI = pars['SBI'][0]
                self.SBI_pars = pars['SBI'][1]
            if 'N_sites' in pars:
                self.N_sites = pars['N_sites'][0]
                self.N_sites_pars = pars['N_sites'][1]
            if 'N_sites_with_DSB' in pars:
                self.N_sites_with_DSB = pars['N_sites_with_DSB'][0]
                self.N_sites_with_DSB_pars = pars['N_sites_with_DSB'][1]
            if 'gamma_par1' in pars:
                self.gamma_par1 = pars['gamma_par1'][0]
                self.gamma_par1_pars = pars['gamma_par1'][1]
            if 'gamma_par2' in pars:
                self.gamma_par2 = pars['gamma_par2'][0]
                self.gamma_par2_pars = pars['gamma_par2'][1]


    @staticmethod
    def linear_func(x, a):
        return a * x

    @staticmethod
    def exp_func(x, a, b):
        return a * (1 - np.exp(-b * x))

    @staticmethod
    def linear_plus_exp_func(x, a, c, d):
        return a * x + c * (1 - np.exp(-d * x))

    @staticmethod
    def quadratic_func(x, a, b, c):
        return a * np.array(x) ** 2 + b * np.array(x) + c

    @staticmethod
    def linquad_func(x, a, b):
        return a * x + b * x ** 2

    @staticmethod
    def survival_func(x, alpha, beta):
        return np.exp(-alpha * x - beta * x ** 2)

    @staticmethod
    def poisson_func(x, a):
        return poisson.pmf(x, a)

    @staticmethod
    def gamma_func(x, a, b):
        rv = gamma(a, b)
        return rv.pdf(x)

    @staticmethod
    def sigmoid_func(x, b, c):
        return 1 / (1 + np.exp(-b * (x - c)))

    @staticmethod
    def get_rsquared(y_pred, y_mean, y):
        SS_res = np.sum((y - y_pred) ** 2)
        SS_tot = np.sum((y - y_mean) ** 2)
        return 1 - (SS_res / SS_tot)

    def fit_compound_poisson(data):
        lambda_poisson = poisson.fit(np.round(data))[0]
        lambda_compound = np.mean(data)
        return lambda_poisson, lambda_compound


